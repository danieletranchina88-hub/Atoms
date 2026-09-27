// Modello molecolare 3D: sfere e bastoncini con legami singoli, doppi e tripli, etichette,
// freccia del momento di dipolo e animazione dei modi normali di vibrazione.

import * as THREE from 'three';
import { textSprite } from './viewer.js';
import { cpkColor, covalentRadius } from '../chem/elementData.js';

const ANG = 1 / 0.52917721090; // Å → bohr

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

/**
 * @param atoms [{Z, xyz}] in bohr
 * @param bonds [{a, b, order}] (order può essere frazionario: 1,5 = aromatico)
 * @param opts { labels: string[] | null, atomColors: string[] | null, scale }
 */
export function buildMolecule(atoms, bonds, opts = {}) {
  const group = new THREE.Group();
  const atomMeshes = [];
  const bondMeshes = [];
  const sphereGeo = new THREE.SphereGeometry(1, 32, 20);
  const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 14, 1, false);
  const light = cssVar('--scene-mode') === 'light';
  const bondColor = new THREE.Color(light ? '#8a93a3' : '#b7bfcc');
  const scale = opts.scale ?? 1;

  atoms.forEach((a, i) => {
    const color = opts.atomColors?.[i] ?? cpkColor(a.Z);
    const mat = new THREE.MeshStandardMaterial({ color, roughness: 0.35, metalness: 0.05 });
    const m = new THREE.Mesh(sphereGeo, mat);
    const r = (0.12 + 0.0034 * covalentRadius(a.Z)) * ANG * scale;
    m.scale.setScalar(r);
    m.position.set(...a.xyz);
    m.userData.radius = r;
    group.add(m);
    atomMeshes.push(m);
  });

  const bondRadius = 0.07 * ANG * scale;
  const matSolid = new THREE.MeshStandardMaterial({ color: bondColor, roughness: 0.5 });
  const matDashed = new THREE.MeshStandardMaterial({ color: bondColor, roughness: 0.5, transparent: true, opacity: 0.45 });

  // direzione perpendicolare al legame, nel piano di un atomo vicino (così i legami doppi stanno nel piano molecolare)
  const perpendicular = (b, pos) => {
    const A = pos[b.a], B = pos[b.b];
    const axis = new THREE.Vector3().subVectors(B, A).normalize();
    let ref = null;
    for (const c of bonds) {
      if (c === b) continue;
      const other = c.a === b.a || c.a === b.b ? c.b : c.b === b.a || c.b === b.b ? c.a : null;
      if (other === null || other === b.a || other === b.b) continue;
      ref = new THREE.Vector3().subVectors(pos[other], A);
      break;
    }
    if (!ref) ref = Math.abs(axis.z) < 0.9 ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(1, 0, 0);
    const perp = ref.sub(axis.clone().multiplyScalar(ref.dot(axis)));
    if (perp.lengthSq() < 1e-8) perp.set(axis.y, -axis.x, 0);
    return perp.normalize();
  };

  const pieces = []; // per l'animazione: { mesh, bond, offset }
  bonds.forEach((b) => {
    const order = b.order;
    let offsets;
    let dashedFrom = Infinity;
    if (order >= 2.75) offsets = [0, 1, -1];
    else if (order >= 2.25) { offsets = [0, 1, -1]; dashedFrom = 2; }        // 2,5
    else if (order >= 1.75) offsets = [0.5, -0.5];
    else if (order >= 1.2) { offsets = [0, 1]; dashedFrom = 1; }               // 1,5 aromatico
    else if (order <= 0.6) { offsets = [0]; dashedFrom = 0; }                  // legame parziale (stato di transizione)
    else offsets = [0];
    const spacing = 0.17 * ANG * scale;
    offsets.forEach((o, k) => {
      const mesh = new THREE.Mesh(cylGeo, k >= dashedFrom ? matDashed : matSolid);
      const r = offsets.length > 1 ? bondRadius * 0.7 : bondRadius;
      mesh.userData = { r };
      group.add(mesh);
      bondMeshes.push(mesh);
      pieces.push({ mesh, b, offset: o * spacing, r });
    });
  });

  const labelSprites = [];
  if (opts.labels) {
    const color = cssVar('--text');
    opts.labels.forEach((txt, i) => {
      if (!txt) return;
      const s = textSprite(txt, color, 30);
      const h = 0.32 * ANG * scale;
      s.scale.set(h * s.userData.aspect, h, 1);
      group.add(s);
      labelSprites.push({ s, i });
    });
  }

  const update = (positions) => {
    const pos = positions.map(p => new THREE.Vector3(...p));
    atomMeshes.forEach((m, i) => m.position.copy(pos[i]));
    for (const pc of pieces) {
      const A = pos[pc.b.a], B = pos[pc.b.b];
      const dir = new THREE.Vector3().subVectors(B, A);
      const len = dir.length();
      const shift = pc.offset ? perpendicular(pc.b, pos).multiplyScalar(pc.offset) : new THREE.Vector3();
      pc.mesh.position.copy(A).add(B).multiplyScalar(0.5).add(shift);
      pc.mesh.scale.set(pc.r, len, pc.r);
      pc.mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.normalize());
    }
    for (const { s, i } of labelSprites) {
      const r = atomMeshes[i].userData.radius;
      s.position.copy(pos[i]).add(new THREE.Vector3(r * 0.9, -r * 0.9, r * 1.2));
    }
  };
  update(atoms.map(a => a.xyz));
  return { group, update, atomMeshes };
}

/** Freccia del momento di dipolo, nella convenzione dei chimici: punta verso l'estremità negativa (δ−). */
export function dipoleArrow(dipoleAU, debye, extent) {
  const v = new THREE.Vector3(...dipoleAU).multiplyScalar(-1); // μ fisico va da − a +; la freccia dei chimici è opposta
  if (v.length() < 1e-6 || debye < 0.02) return null;
  const len = Math.min(extent * 1.1, (1.2 + 0.45 * debye) * ANG);
  const dir = v.normalize();
  const color = cssVar('--accent');
  const group = new THREE.Group();
  const origin = dir.clone().multiplyScalar(-len / 2);
  const arrow = new THREE.ArrowHelper(dir, origin, len, color, 0.22 * len, 0.1 * len);
  // sempre visibile, anche dentro gli atomi
  arrow.traverse(o => { if (o.material) { o.material.depthTest = false; o.material.transparent = true; o.renderOrder = 20; } });
  group.add(arrow);
  // barretta trasversale della "croce" all'estremità positiva
  const perp = Math.abs(dir.z) < 0.9 ? new THREE.Vector3(0, 0, 1).cross(dir).normalize() : new THREE.Vector3(1, 0, 0).cross(dir).normalize();
  const c = origin.clone().add(dir.clone().multiplyScalar(len * 0.12));
  const geo = new THREE.BufferGeometry().setFromPoints([c.clone().add(perp.clone().multiplyScalar(len * 0.08)), c.clone().add(perp.clone().multiplyScalar(-len * 0.08))]);
  const bar = new THREE.Line(geo, new THREE.LineBasicMaterial({ color, depthTest: false, transparent: true }));
  bar.renderOrder = 20;
  group.add(bar);
  return group;
}

/** Colore per il potenziale elettrostatico: rosso (negativo) → bianco → blu (positivo). */
export function espColor(v, vmax, out, offset) {
  const t = Math.max(-1, Math.min(1, v / vmax));
  let r, g, b;
  if (t < 0) { r = 1; g = 1 + t * 0.85; b = 1 + t * 0.95; } else { r = 1 - t * 0.95; g = 1 - t * 0.6; b = 1; }
  out[offset] = r; out[offset + 1] = g; out[offset + 2] = b;
}

/** Colore per una carica parziale: rosso (negativa) → grigio → blu (positiva). */
export function chargeColor(q) {
  const t = Math.max(-1, Math.min(1, q / 0.8));
  const c = new THREE.Color();
  if (t < 0) c.setRGB(0.85, 0.85 + 0.75 * t, 0.85 + 0.8 * t);
  else c.setRGB(0.85 - 0.75 * t, 0.85 - 0.45 * t, 0.9);
  return `#${c.getHexString()}`;
}
