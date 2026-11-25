"""
Exporta os resultados de fragmentação e conectividade (Artigo 6, BHSL) para o visualizador 3D em Three.js.
Mesmo esquema do visualizador do ICEA (Artigo 3), inspirado em linaq11/lidar-mangrove-canopy.

Entradas (data/; geradas por fragconn_analise.py, ver o pacote Zenodo 10.5281/zenodo.24137643):
  A_mapas.npz        m14, m24 (0 matriz, 1 núcleo, 2 borda, 3 ilhota), lost, gain, inside
  C_nos_2024.npz     rótulos dos nós 2024 (fragmentos >= 10 ha) e transform (EPSG:31981, 30 m)
  C01_dPC_2024_NAT_1500m.csv, C02_importancia_2014_perda.csv, C03_top_conectores.csv
  D01_processos_celulas_5km.csv, D_grid.npz
  A01_metricas_anuais.csv, B01_eca_anual.csv, E01_taxas_perda.csv, numeros_texto.json, parametros.json
  CESL_V01_limite_corredor.shp (limite do corredor, EPSG:4326)

Saída:
  docs/index.html (arquivo único; publicado pelo GitHub Pages)

Uso (na raiz do repositório): python scripts/exportar_visualizacao_3d_fragconn.py [--resultados DIR] [--limite SHP] [--saida HTML]
"""
import argparse
import base64
import json
import math
import os
import sys
import zlib

import numpy as np
import pandas as pd
import shapefile
from scipy import ndimage as ndi

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CORR = os.path.join(BASE, "data", "CESL_V01_limite_corredor.shp")

ap = argparse.ArgumentParser(description="Gera o visualizador 3D de fragmentação e conectividade (HTML único).")
ap.add_argument("--resultados", default=os.path.join(BASE, "data"))
ap.add_argument("--limite", default=_CORR, help="shapefile do limite do corredor (EPSG:4326)")
ap.add_argument("--modelo", default=os.path.join(BASE, "web", "visualizador_fragconn_3d.template.html"))
ap.add_argument("--saida", default=os.path.join(BASE, "docs", "index.html"))
ARGS = ap.parse_args()
RES = ARGS.resultados

PIX = 30.0
BLOCO = 8            # visão geral: 8 x 8 células (240 m)
JANELA = 192         # janelas de detalhe: 192 x 192 células (5,76 km)
D_CAP = 1000.0       # distância à borda saturada em 1 km (última faixa de borda do artigo)
D_PASSO = 4.0        # codificação: 4 m por nível (0..250)
S8 = np.ones((3, 3), bool)
CODES = {'agregacao': 1, 'atricao': 2, 'criacao': 3, 'dissecacao': 4, 'ampliacao': 5, 'fragmentacao': 6,
         'perfuracao': 7, 'encolhimento': 8, 'estavel': 9}


def b64z(a):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(a, dtype=np.uint8).tobytes(), 9)).decode("ascii")


def enc_log(v):        # importância (%) em escala log: 1e-4 .. 100 -> 1..250; 0 = sem nó
    v = np.asarray(v, float)
    e = np.where(v > 0, np.clip(np.rint((np.log10(np.maximum(v, 1e-4)) + 4) / 6 * 249) + 1, 1, 250), 0)
    return e.astype(np.uint8)


def log(*a):
    print(*a, flush=True)


def main():
    mp = np.load(os.path.join(RES, "A_mapas.npz"))
    m = {2014: mp["m14"], 2024: mp["m24"]}
    lost, gain, inside = mp["lost"], mp["gain"], mp["inside"]
    H, W = inside.shape
    cn = np.load(os.path.join(RES, "C_nos_2024.npz"))
    tr = cn["transform"]            # a, b, c, d, e, f (affine)
    x0, y0 = float(tr[2]), float(tr[5])
    log(f"grade {W} x {H}, origem {x0:.0f} {y0:.0f}")

    # ---------------- distância à borda e nós por ano ----------------
    dist, nodes, info = {}, {}, {}
    c01 = pd.read_csv(os.path.join(RES, "C01_dPC_2024_NAT_1500m.csv")).set_index("node")
    c02 = pd.read_csv(os.path.join(RES, "C02_importancia_2014_perda.csv")).set_index("node14")
    c03 = pd.read_csv(os.path.join(RES, "C03_top_conectores.csv"))
    for ano in (2014, 2024):
        hab = m[ano] > 0
        matriz = inside & ~hab
        dist[ano] = np.minimum(ndi.distance_transform_edt(~matriz, sampling=PIX), D_CAP).astype(np.float32)
        dist[ano][~hab] = 0
        # nós = fragmentos >= 10 ha, mesma rotulagem do fragconn_analise.py (vizinhança 8)
        lab, n = ndi.label(hab, structure=S8)
        areas = np.bincount(lab.ravel())[1:] * PIX * PIX / 1e4
        keep = np.where(areas >= 10.0)[0] + 1
        remap = np.zeros(n + 1, np.int32); remap[keep] = np.arange(1, len(keep) + 1)
        nd = remap[lab]; del lab
        if ano == 2024:
            assert (nd == cn["nodes"]).all(), "rotulagem dos nós 2024 difere da do C_nos_2024.npz"
        nodes[ano] = nd
        a = areas[keep - 1]
        # ligações do grafo planar (mesma regra do script de análise)
        D, (iy, ix) = ndi.distance_transform_edt(nd == 0, sampling=PIX, return_indices=True)
        A = nd[iy, ix]; A[~inside] = 0; del iy, ix
        pa, pb, pdd = [], [], []
        for s1, s2 in (((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
                       ((slice(None, -1), slice(None)), (slice(1, None), slice(None)))):
            a1, a2 = A[s1], A[s2]
            k = (a1 != a2) & (a1 > 0) & (a2 > 0)
            pa.append(a1[k]); pb.append(a2[k]); pdd.append((D[s1][k] + D[s2][k]).astype(np.float32))
        del A, D
        df = pd.DataFrame({"i": np.minimum(np.concatenate(pa), np.concatenate(pb)),
                           "j": np.maximum(np.concatenate(pa), np.concatenate(pb)),
                           "d": np.concatenate(pdd)}).groupby(["i", "j"], as_index=False)["d"].min()
        df["d"] = df["d"].clip(lower=PIX)
        idx = np.arange(1, len(a) + 1)
        cy, cx = np.array(ndi.center_of_mass(np.ones_like(nd, np.uint8), nd, idx)).T
        if ano == 2024:
            imp = c01.loc[idx, "dPC"].values
            con = c01.loc[idx, "dPCconnector"].values
            intra = c01.loc[idx, "dPCintra"].values
            flux = c01.loc[idx, "dPCflux"].values
        else:
            imp = c02.loc[idx, "imp_intra_flux_pct"].values
            con = np.full(len(a), np.nan); intra = flux = con
        rot = {int(r.node): r.rotulo for r in c03.itertuples()} if ano == 2024 else {}
        perda = c02.loc[idx, "frac_lost_by_2024"].values if ano == 2014 else np.full(len(a), np.nan)
        info[ano] = {
            "nos": [[round(float(cx[k]), 1), round(float(cy[k]), 1), round(float(a[k]), 2), float(f"{imp[k]:.4g}"),
                     None if np.isnan(con[k]) else float(f"{con[k]:.4g}"),
                     None if np.isnan(intra[k]) else float(f"{intra[k]:.4g}"),
                     None if np.isnan(flux[k]) else float(f"{flux[k]:.4g}"),
                     None if np.isnan(perda[k]) else round(float(perda[k]), 3),
                     rot.get(k + 1, "")] for k in range(len(a))],
            "lig": df[["i", "j", "d"]].astype({"i": int, "j": int}).assign(d=lambda t: t.d.round(0).astype(int)).values.tolist(),
        }
        log(f"{ano}: {len(a)} nós, {len(df)} ligações")

    # ---------------- processos de Bogaert por célula de 5 km ----------------
    grid = np.load(os.path.join(RES, "D_grid.npz"))["grid"]
    d01 = pd.read_csv(os.path.join(RES, "D01_processos_celulas_5km.csv"))
    lut = np.zeros(grid.max() + 1, np.uint8)
    for r in d01.itertuples():
        lut[int(r.cell)] = CODES[r.processo]
    proc = lut[grid]

    # ---------------- visão geral (blocos) ----------------
    Hb, Wb = H // BLOCO, W // BLOCO
    def bl(a):
        return a[: Hb * BLOCO, : Wb * BLOCO].reshape(Hb, BLOCO, Wb, BLOCO)
    nin = bl(inside).sum(axis=(1, 3)).astype(np.float32)
    den = np.maximum(nin, 1)
    def frac(mask):
        return np.clip(np.rint(bl(mask & inside).sum(axis=(1, 3)) / den * 250), 0, 250).astype(np.uint8)
    geral = {"dentro": b64z(np.clip(np.rint(nin / (BLOCO * BLOCO) * 250), 0, 250)),
             "perda": b64z(frac(lost)), "ganho": b64z(frac(gain))}
    pc = proc[BLOCO // 2:: BLOCO, BLOCO // 2:: BLOCO][:Hb, :Wb]
    geral["proc"] = b64z(pc)
    for ano in (2014, 2024):
        prof = bl(dist[ano]).sum(axis=(1, 3)) / den
        geral[f"prof{ano}"] = b64z(np.clip(np.rint(prof / D_PASSO), 0, 250))
        for c, nome in ((1, "nucleo"), (2, "borda"), (3, "ilhota")):
            geral[f"{nome}{ano}"] = b64z(frac(m[ano] == c))
        impv = np.zeros(len(info[ano]["nos"]) + 1); impv[1:] = [n_[3] for n_ in info[ano]["nos"]]
        impmax = bl(impv[nodes[ano]]).max(axis=(1, 3))
        geral[f"imp{ano}"] = b64z(enc_log(impmax))

    # ---------------- janelas de detalhe ----------------
    def janela_centro(r, c):
        r0 = int(np.clip(r - JANELA // 2, 0, H - JANELA)); c0 = int(np.clip(c - JANELA // 2, 0, W - JANELA))
        return r0, c0
    escolhas = []
    # 1. conector C3 (fragmento pequeno de alta importância como conector)
    c3 = c03.sort_values("rotulo").iloc[2]
    no = int(c3.node); n24 = info[2024]["nos"][no - 1]
    escolhas.append(("Conector C3: fragmento de %s ha" % f"{c3.area_ha:.0f}".replace(".", ","), *janela_centro(n24[1], n24[0]), "conector"))
    # 2. célula de dissecação com maior queda de meff
    dd = d01[(d01.processo == "dissecacao") & (d01.a14 > 300)].sort_values("d_meff_pct").iloc[0]
    escolhas.append(("Dissecação (célula de 5 km)", *janela_centro(dd.row + 83, dd.col + 83), "dissecacao"))
    # 3. maior perda 2014–2024 e 4. núcleo mais estável
    melhor_perda, melhor_nucleo = (-1, 0, 0), (-1, 0, 0)
    cs_l = np.pad(lost.astype(np.int32).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    cs_c = np.pad((m[2024] == 1).astype(np.int32).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    cs_i = np.pad(inside.astype(np.int32).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    soma = lambda cs, r, c: cs[r + JANELA, c + JANELA] - cs[r, c + JANELA] - cs[r + JANELA, c] + cs[r, c]
    for r in range(0, H - JANELA, 48):
        for c in range(0, W - JANELA, 48):
            if soma(cs_i, r, c) < 0.95 * JANELA * JANELA:
                continue
            vl, vc = soma(cs_l, r, c), soma(cs_c, r, c)
            if vl > melhor_perda[0]: melhor_perda = (vl, r, c)
            if vc > melhor_nucleo[0]: melhor_nucleo = (vc, r, c)
    escolhas.append(("Maior perda de vegetação 2014–2024", melhor_perda[1], melhor_perda[2], "perda"))
    escolhas.append(("Núcleo extenso preservado", melhor_nucleo[1], melhor_nucleo[2], "nucleo"))

    janelas = []
    for nome, r0, c0, crit in escolhas:
        sl = (slice(r0, r0 + JANELA), slice(c0, c0 + JANELA))
        ins = inside[sl]
        jd = {"nome": nome, "criterio": crit, "r0": int(r0), "c0": int(c0), "n": JANELA, "dados": {}}
        mud = np.zeros(ins.shape, np.uint8); mud[lost[sl]] = 1; mud[gain[sl]] = 2; mud[~ins] = 255
        jd["mud"] = b64z(mud)
        jd["proc"] = int(proc[r0 + JANELA // 2, c0 + JANELA // 2])
        for ano in (2014, 2024):
            mo = m[ano][sl].astype(np.uint8).copy(); mo[~ins] = 255
            impv = np.zeros(len(info[ano]["nos"]) + 1); impv[1:] = [n_[3] for n_ in info[ano]["nos"]]
            jd["dados"][str(ano)] = {
                "morf": b64z(mo),
                "prof": b64z(np.clip(np.rint(dist[ano][sl] / D_PASSO), 0, 250)),
                "imp": b64z(enc_log(impv[nodes[ano][sl]])),
                "no": b64z(np.minimum(nodes[ano][sl], 65535).astype("<u2").view(np.uint8)),
            }
        jd["perda_ha"] = round(float(lost[sl].sum() * PIX * PIX / 1e4), 1)
        janelas.append(jd)
        log(f"  janela: {nome} | lin {r0} col {c0} | perda {jd['perda_ha']} ha")

    # ---------------- limite do corredor (4326 -> 31981) ----------------
    limite = []
    if ARGS.limite and os.path.exists(ARGS.limite):
        from pyproj import Transformer
        tf = Transformer.from_crs("EPSG:4326", "EPSG:31981", always_xy=True)
        for sh in shapefile.Reader(ARGS.limite).shapes():
            partes = list(sh.parts) + [len(sh.points)]
            for a_, b_ in zip(partes[:-1], partes[1:]):
                xs, ys = tf.transform(*zip(*sh.points[a_:b_]))
                limite.append([[round((x - x0) / PIX, 1), round((y0 - y) / PIX, 1)] for x, y in zip(xs, ys)])
    else:
        log("[aviso] limite do corredor não encontrado; a página sai sem o contorno")

    # ---------------- números do painel (lidos dos CSV/JSON) ----------------
    a01 = pd.read_csv(os.path.join(RES, "A01_metricas_anuais.csv"))
    b01 = pd.read_csv(os.path.join(RES, "B01_eca_anual.csv"))
    e01 = pd.read_csv(os.path.join(RES, "E01_taxas_perda.csv"))
    nt = json.load(open(os.path.join(RES, "numeros_texto.json"), encoding="utf-8"))
    nat = a01[a01.habitat == "NAT"].set_index("year")
    eca = b01[(b01.habitat == "NAT")].pivot(index="year", columns="dist_m", values="ECA_ha")
    numeros = {
        "serie": {"anos": nat.index.tolist(), "habitat_ha": nat.area_ha.round(1).tolist(),
                  "eca1500_ha": eca[1500].round(1).tolist()},
        "habitat": {str(y): round(float(nat.loc[y, "area_ha"]), 1) for y in (2014, 2024)},
        "nucleo_pct": {str(y): round(float(nat.loc[y, "core_pct_hab"]), 1) for y in (2014, 2024)},
        "nos": {str(y): int(nat.loc[y, "NP_ge10ha"]) for y in (2014, 2024)},
        "eca1500": {str(y): round(float(eca.loc[y, 1500]), 1) for y in (2014, 2024)},
        "eca_var_pct": {str(d): round(100 * (eca.loc[2024, d] - eca.loc[2014, d]) / eca.loc[2014, d], 1) for d in eca.columns},
        "perda_bruta_ha": round(nt["perda_bruta_ha"], 1),
        "savana_pct_perda": round(nt["perda_bruta_savana_pct_do_total"], 1),
        "perda_borda": e01[e01.variavel == "dist_c"][["classe", "loss_rate_pct"]].round(2).values.tolist(),
        "processos": nt["processos"], "celulas": nt["celulas"],
    }

    dados = {
        "meta": {"W": W, "H": H, "pix": PIX, "bloco": BLOCO, "Wb": Wb, "Hb": Hb, "d_cap": D_CAP, "d_passo": D_PASSO,
                 "d_focal": 1500, "x0": x0, "y0": y0, "crs": "EPSG:31981"},
        "geral": geral, "rede": info, "janelas": janelas, "limite": limite, "numeros": numeros,
    }
    js = json.dumps(dados, ensure_ascii=False, separators=(",", ":"))
    html = open(ARGS.modelo, encoding="utf-8").read().replace("/*__DADOS__*/null", js)
    os.makedirs(os.path.dirname(ARGS.saida), exist_ok=True)
    open(ARGS.saida, "w", encoding="utf-8", newline="\n").write(html)
    log(f"visão geral {Wb}x{Hb} | {ARGS.saida} ({os.path.getsize(ARGS.saida) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
