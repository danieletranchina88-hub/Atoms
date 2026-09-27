// Simmetria molecolare: elementi di simmetria, gruppo puntuale (notazione di Schoenflies)
// e numero di simmetria rotazionale σ (necessario per la funzione di partizione rotazionale).

const TOL = 0.08; // bohr

function centerOfMass(atoms, masses) {
  const c = [0, 0, 0];
  let M = 0;
  atoms.forEach((a, i) => { for (let k = 0; k < 3; k++) c[k] += masses[i] * a.xyz[k]; M += masses[i]; });
  return c.map(v => v / M);
}

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => { const l = norm(a); return l > 1e-9 ? [a[0] / l, a[1] / l, a[2] / l] : null; };

/** Matrice di rotazione di angolo θ attorno all'asse unitario u (formula di Rodrigues). */
function rotation(u, th) {
  const c = Math.cos(th), s = Math.sin(th), C = 1 - c;
  const [x, y, z] = u;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}
const reflection = (u) => [
  [1 - 2 * u[0] * u[0], -2 * u[0] * u[1], -2 * u[0] * u[2]],
  [-2 * u[1] * u[0], 1 - 2 * u[1] * u[1], -2 * u[1] * u[2]],
  [-2 * u[2] * u[0], -2 * u[2] * u[1], 1 - 2 * u[2] * u[2]],
];
const mul = (A, B) => A.map((r, i) => [0, 1, 2].map(j => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const apply = (M, v) => [dot(M[0], v), dot(M[1], v), dot(M[2], v)];

/** Verifica che l'operazione M trasformi la molecola in sé stessa. */
function isSymmetry(pts, Zs, M) {
  for (let i = 0; i < pts.length; i++) {
    const q = apply(M, pts[i]);
    let found = false;
    for (let j = 0; j < pts.length; j++) {
      if (Zs[j] !== Zs[i]) continue;
      if (Math.hypot(q[0] - pts[j][0], q[1] - pts[j][1], q[2] - pts[j][2]) < TOL) { found = true; break; }
    }
    if (!found) return false;
  }
  return true;
}

/** Tensore d'inerzia e assi principali. */
export function inertia(atoms, masses) {
  const c = centerOfMass(atoms, masses);
  const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  atoms.forEach((a, i) => {
    const r = [a.xyz[0] - c[0], a.xyz[1] - c[1], a.xyz[2] - c[2]];
    const r2 = dot(r, r);
    for (let p = 0; p < 3; p++) for (let q = 0; q < 3; q++) I[p][q] += masses[i] * ((p === q ? r2 : 0) - r[p] * r[q]);
  });
  // diagonalizzazione di Jacobi 3×3
  const A = I.map(r => [...r]);
  const V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) off += A[p][q] ** 2;
    if (off < 1e-24) break;
    for (let p = 0; p < 3; p++) {
      for (let q = p + 1; q < 3; q++) {
        if (Math.abs(A[p][q]) < 1e-30) continue;
        const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
        const c2 = Math.cos(th), s2 = Math.sin(th);
        for (let k = 0; k < 3; k++) {
          const akp = A[k][p], akq = A[k][q];
          A[k][p] = c2 * akp - s2 * akq;
          A[k][q] = s2 * akp + c2 * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = A[p][k], aqk = A[q][k];
          A[p][k] = c2 * apk - s2 * aqk;
          A[q][k] = s2 * apk + c2 * aqk;
        }
        for (let k = 0; k < 3; k++) {
          const vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c2 * vkp - s2 * vkq;
          V[k][q] = s2 * vkp + c2 * vkq;
        }
      }
    }
  }
  const moments = [A[0][0], A[1][1], A[2][2]];
  const axes = [0, 1, 2].map(k => [V[0][k], V[1][k], V[2][k]]);
  const order = [0, 1, 2].sort((a, b) => moments[a] - moments[b]);
  return { center: c, moments: order.map(k => moments[k]), axes: order.map(k => axes[k]) };
}

/**
 * Analisi di simmetria completa.
 * Restituisce { pointGroup, sigma (numero di simmetria rotazionale), linear, elements }.
 */
export function analyzeSymmetry(atoms, masses) {
  const { center, moments, axes } = inertia(atoms, masses);
  const pts = atoms.map(a => [a.xyz[0] - center[0], a.xyz[1] - center[1], a.xyz[2] - center[2]]);
  const Zs = atoms.map(a => a.Z);
  const N = atoms.length;
  if (N === 1) return { pointGroup: 'K_h', sigma: 1, linear: false, atom: true, elements: {} };
  const linear = moments[0] < 1e-3 * Math.max(moments[2], 1e-9);
  const inversion = isSymmetry(pts, Zs, [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);
  if (linear) {
    return {
      pointGroup: inversion ? 'D∞h' : 'C∞v', sigma: inversion ? 2 : 1, linear: true,
      elements: { inversion },
    };
  }
  // assi candidati: assi principali, direzioni degli atomi, punti medi, normali di coppie di atomi
  const cand = [...axes];
  for (let i = 0; i < N; i++) {
    const u = unit(pts[i]);
    if (u) cand.push(u);
    for (let j = i + 1; j < N; j++) {
      if (Zs[i] !== Zs[j]) continue;
      const m = unit([(pts[i][0] + pts[j][0]) / 2, (pts[i][1] + pts[j][1]) / 2, (pts[i][2] + pts[j][2]) / 2]);
      if (m) cand.push(m);
      const c = unit(cross(pts[i], pts[j]));
      if (c) cand.push(c);
      const d = unit([pts[i][0] - pts[j][0], pts[i][1] - pts[j][1], pts[i][2] - pts[j][2]]);
      if (d) cand.push(d);
    }
  }
  // elimina duplicati (anche con verso opposto)
  const uniqueAxes = [];
  for (const u of cand) {
    if (!uniqueAxes.some(v => Math.abs(Math.abs(dot(u, v)) - 1) < 1e-4)) uniqueAxes.push(u);
  }
  // rotazioni proprie Cn (n = 2…8) e piani di simmetria
  const rotations = [];
  const cnAxes = [];
  const planes = [];
  for (const u of uniqueAxes) {
    let maxN = 1;
    for (let n = 8; n >= 2; n--) {
      if (isSymmetry(pts, Zs, rotation(u, 2 * Math.PI / n))) { maxN = n; break; }
    }
    if (maxN > 1) {
      cnAxes.push({ axis: u, n: maxN });
      for (let k = 1; k < maxN; k++) rotations.push(rotation(u, 2 * Math.PI * k / maxN));
    }
    if (isSymmetry(pts, Zs, reflection(u))) planes.push(u);
  }
  // numero di simmetria = ordine del sottogruppo delle rotazioni proprie (identità inclusa)
  const distinct = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]];
  for (const R of rotations) {
    if (!distinct.some(D => D.every((row, i) => row.every((v, j) => Math.abs(v - R[i][j]) < 1e-4)))) distinct.push(R);
  }
  const sigma = distinct.length;
  // asse principale e gruppo puntuale (diagramma di flusso standard)
  const nMax = cnAxes.reduce((m, a) => Math.max(m, a.n), 1);
  const c2count = cnAxes.filter(a => a.n % 2 === 0).length;
  const highOrder = cnAxes.filter(a => a.n >= 3).length;
  let pg;
  if (highOrder >= 2) {
    if (cnAxes.some(a => a.n === 5)) pg = inversion ? 'Ih' : 'I';
    else if (cnAxes.some(a => a.n === 4)) pg = inversion ? 'Oh' : 'O';
    else pg = planes.length ? (inversion ? 'Th' : 'Td') : 'T';
  } else if (nMax === 1) {
    pg = planes.length ? 'Cs' : inversion ? 'Ci' : 'C1';
  } else {
    const main = cnAxes.find(a => a.n === nMax).axis;
    const perpC2 = cnAxes.filter(a => a.n % 2 === 0 && Math.abs(dot(a.axis, main)) < 1e-3).length;
    const sigmaH = planes.some(p => Math.abs(Math.abs(dot(p, main)) - 1) < 1e-3);
    const sigmaV = planes.filter(p => Math.abs(dot(p, main)) < 1e-3).length;
    if (perpC2 >= 1) {
      pg = sigmaH ? `D${nMax}h` : sigmaV ? `D${nMax}d` : `D${nMax}`;
    } else if (sigmaH) {
      pg = `C${nMax}h`;
    } else if (sigmaV) {
      pg = `C${nMax}v`;
    } else {
      // S2n
      const s2n = isSymmetry(pts, Zs, mul(reflection(main), rotation(main, Math.PI / nMax)));
      pg = s2n ? `S${2 * nMax}` : `C${nMax}`;
    }
  }
  return {
    pointGroup: pg, sigma, linear: false,
    elements: { inversion, planes: planes.length, axes: cnAxes.map(a => a.n), c2: c2count },
    principal: { moments, axes, center },
  };
}
