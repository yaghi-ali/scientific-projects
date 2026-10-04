"""TD4 - EDFA. Python 3 et NumPy. Exécution : python correction_edfa.py.
Produit courbes_edfa.html dans le dossier du script.
Les formules ASE et NF sont celles imposées par le TD, pour G > 1.
"""
from pathlib import Path
import numpy as np

h, c = 6.62607015e-34, 299792458.0
lambda_s, lambda_p = 1550e-9, 980e-9
nu_s, nu_p = c / lambda_s, c / lambda_p
sigma_es, sigma_as = 3e-25, 1.5e-25
sigma_ap, sigma_ep = 2e-25, 0.0
tau, Ntot, Gamma, Aeff = 10e-3, 1e25, 0.8, 50e-12
alpha = 0.005 * np.log(10) / 10  # dB/m -> coefficient de puissance en m^-1
nsp, delta_lambda = 1.2, 0.1e-9
delta_nu = c * delta_lambda / lambda_s**2
C_ASE = 2 * nsp * h * nu_s * delta_nu

def simulate(Ps0_dBm=-30, Pp0=62e-3, L=18.0, dz=1e-4):
    if min(Pp0, L, dz) <= 0 or not np.isfinite([Ps0_dBm,Pp0,L,dz]).all():
        raise ValueError('Finite parameters and positive pump, length and step required.')
    Nz = int(np.ceil(L / dz))
    z = np.linspace(0, L, Nz + 1)
    dz = z[1] - z[0]
    Ps, Pp = np.empty(Nz + 1), np.empty(Nz + 1)
    Ps[0], Pp[0] = 1e-3 * 10**(Ps0_dBm / 10), Pp0
    for i in range(Nz):
        phi_s = Gamma * Ps[i] / (h * nu_s * Aeff)
        phi_p = Gamma * Pp[i] / (h * nu_p * Aeff)
        N2 = Ntot * (sigma_ap * phi_p + sigma_as * phi_s) / (
            (sigma_ap + sigma_ep) * phi_p
            + (sigma_as + sigma_es) * phi_s + 1 / tau)
        N1 = Ntot - N2
        gs = Gamma * (sigma_es * N2 - sigma_as * N1) - alpha
        ap = Gamma * (sigma_ap * N1 - sigma_ep * N2) + alpha
        Ps[i+1] = Ps[i] * (1 + gs * dz)
        Pp[i+1] = Pp[i] * (1 - ap * dz)
        if Ps[i+1] <= 0 or Pp[i+1] <= 0:
            raise ValueError('Pas Euler trop grand : diminuer dz.')
    G = Ps / Ps[0]
    valid = G > 1
    ASE = np.full_like(G, np.nan)
    OSNR = np.full_like(G, np.nan)
    NF = np.full_like(G, np.nan)
    ASE[valid] = C_ASE * (G[valid] - 1)
    OSNR[valid] = 10 * np.log10(Ps[valid] / ASE[valid])
    NF[valid] = 10 * np.log10(2 * nsp * G[valid] / (G[valid] - 1))
    return dict(z=z, Ps=Ps, Pp=Pp, G=G, GdB=10*np.log10(G),
                ASE=ASE, OSNR=OSNR, NF=NF, k=int(np.argmax(G)))

# SVG : courbes consultables sans installer une bibliothèque graphique.
def plot(title, ylabel, series, ymin, ymax):
    W, H, left, top, width, height = 860, 340, 72, 36, 755, 245
    parts = [f'<section><h2>{title}</h2><svg viewBox="0 0 {W} {H}" role="img" aria-label="{title}">',
             f'<text x="8" y="20">{ylabel}</text>']
    for value in np.linspace(ymin, ymax, 6):
        y = top + height * (ymax-value)/(ymax-ymin)
        parts.append(f'<line x1="{left}" y1="{y}" x2="{left+width}" y2="{y}" stroke="#ddd"/><text x="{left-8}" y="{y+4}" text-anchor="end">{value:.1f}</text>')
    for value in range(0, 19, 3):
        x = left + width*value/18
        parts.append(f'<text x="{x}" y="{top+height+22}" text-anchor="middle">{value}</text>')
    parts.append(f'<text x="{left+width/2}" y="{H-5}" text-anchor="middle">Longueur z (m)</text>')
    for n, (label, z, values, color) in enumerate(series):
        # Décimation uniquement pour l'affichage ; calculs sur tous les points.
        indices = np.unique(np.r_[np.arange(0,len(z),max(1,len(z)//2400)),len(z)-1])
        points = []
        for i in indices:
            v = values[i]
            if np.isfinite(v) and ymin <= v <= ymax:
                points.append(f'{left+width*z[i]/18:.2f},{top+height*(ymax-v)/(ymax-ymin):.2f}')
            elif points:
                parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="'+ ' '.join(points)+'"/>')
                points=[]
        if points:
            parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="'+ ' '.join(points)+'"/>')
        parts.append(f'<text x="{left+15+n*230}" y="{top+16}" fill="{color}">{label}</text>')
    return ''.join(parts)+'</svg></section>'

def main():
    cases = [('Cas initial',-30,.062),('Préampli',-40,.050),('Booster',-5,.050)]
    results = [(name,simulate(dbm,pump)) for name,dbm,pump in cases]
    html = ['<!doctype html><html lang="fr"><meta charset="utf-8"><title>TD4 EDFA - Courbes et résultats</title><style>body{font:16px system-ui;max-width:1000px;margin:35px auto;padding:0 20px;color:#172b45;background:#f4f7fb}section{background:white;padding:20px;margin:20px 0;border-radius:12px}svg{width:100%;font:13px system-ui}table{border-collapse:collapse;width:100%}td,th{padding:12px;border-bottom:1px solid #ddd;text-align:right}h1,h2{color:#173c69}p{line-height:1.6}</style><h1>TD4 — Simulation d’un EDFA</h1><p>Euler explicite, L = 18 m, Δz = 0,0001 m. Puissances en watts dans les calculs. Bande OSNR : 0,1 nm. ASE et NF selon les expressions du TD, uniquement pour G &gt; 1.</p><section><h2>Résultats à la longueur optimale</h2><table><tr><th>Cas</th><th>L opt. (m)</th><th>Gain (dB)</th><th>Signal (mW)</th><th>OSNR (dB)</th><th>NF (dB)</th></tr>']
    for name,r in results:
        k=r['k']
        vals=(r['z'][k],r['GdB'][k],r['Ps'][k]*1e3,r['OSNR'][k],r['NF'][k])
        html.append('<tr><td>'+name+'</td>'+''.join(f'<td>{v:.4f}</td>' for v in vals)+'</tr>')
        print(name, 'Lopt,Gmax,Ps_mW,OSNR,NF=', vals, 'gain_18m=',r['GdB'][-1])
    html.append('</table></section>')
    r=results[0][1]
    html.append(plot('1. Cas initial : signal et pompe','Puissance (dBm)',[
        ('Signal',r['z'],10*np.log10(r['Ps']/1e-3),'#006bb6'),
        ('Pompe',r['z'],10*np.log10(r['Pp']/1e-3),'#df6b19')],-70,20))
    html.append(plot('2. Cas initial : gain cumulatif','Gain (dB)',[('Gain',r['z'],r['GdB'],'#006bb6')],0,40))
    html.append(plot('3. Préampli et booster : comparaison des gains','Gain (dB)',[
        (name,q['z'],q['GdB'],color) for (name,q),color in zip(results[1:],['#006bb6','#df6b19'])],-20,40))
    q=results[1][1]
    html.append(plot('4. Préampli : OSNR dans une bande de 0,1 nm','OSNR (dB)',[('OSNR (G > 1)',q['z'],q['OSNR'],'#087f70')],14,45))
    html.append('<section><h2>Interprétation</h2><p>Le gain augmente tant que la pompe entretient une inversion suffisante. Après le maximum, la fibre réabsorbe le signal. Le booster sature davantage le gain et atteint son maximum sur une fibre plus courte.</p><p>Dans la formule imposée, OSNR = [Ps(0)/(2 nsp h νs Δν)] × G/(G−1). Elle tend vers 14,16 dB à fort gain pour le préampli. La remontée lorsque G redescend vers 1 est une limite du modèle global : il ne décrit pas correctement une section absorbante placée après une section amplificatrice. Les points G ≤ 1 sont exclus. Les divergences au voisinage de G = 1 peuvent dépasser le cadre du graphique.</p><p>La NF utilise strictement F = 2 nsp G/(G−1), fourni par le TD. Sa limite à fort gain vaut 3,802 dB. L’ASE est calculée après propagation et ne participe pas à la saturation dans ce modèle.</p></section></html>')
    path=Path(__file__).with_name('courbes_edfa.html')
    path.write_text(''.join(html),encoding='utf-8')
    print('Courbes :',path)

if __name__ == '__main__':
    main()
