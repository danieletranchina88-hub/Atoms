// Kit molecolare. Si costruisce come la scatola di palline.
// La geometria è quella che metti tu; orbitali, ibridazione e densità
// vengono dal calcolo Hartree–Fock (6-31G), non da un disegno.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { request, cancelAll } from '../chem/moleculeClient.js';
import { analyzeStructure } from '../chem/structure.js';
import { cpkColor, covalentRadius } from '../chem/elementData.js';
import { marchingCubes } from '../render/marching.js';
import { ELEMENTS } from '../physics/elements.js';

const BOHR = 1.8897259886; // Å → bohr
const PALETTE = [
  { Z: 1, name: 'H' }, { Z: 6, name: 'C' }, { Z: 7, name: 'N' }, { Z: 8, name: 'O' },
  { Z: 9, name: 'F' }, { Z: 15, name: 'P' }, { Z: 16, name: 'S' }, { Z: 17, name: 'Cl' },
  { Z: 35, name: 'Br' }, { Z: 53, name: 'I' },
];

function snapDistance(za, zb) {
  return 1.35 * (covalentRadius(za) + covalentRadius(zb)) / 100;
}

export function mountKit(root) {
  root.innerHTML = `
    <div class="kit">
      <aside class="kit-box">
        <h2>Scatola</h2>
        <p class="kit-note">Scegli un atomo e clicca sul piano per posarlo. Avvicinali: il legame scatta se la distanza è entro 1,35 volte la somma dei raggi covalenti di Pyykkö.</p>
        <div class="kit-atoms" id="kit-atoms"></div>
        <div class="kit-tools">
          <button type="button" id="kit-bond">Ordine legame</button>
          <button type="button" id="kit-del">Elimina</button>
          <button type="button" id="kit-clear">Svuota</button>
        </div>
        <button type="button" class="kit-calc" id="kit-calc">Calcola Hartree–Fock</button>
        <div id="kit-status" class="kit-note"></div>
        <div id="kit-info"></div>
      </aside>
      <div class="kit-stage" id="kit-stage"></div>
    </div>`;

  const atomsBox = root.querySelector('#kit-atoms');
  atomsBox.innerHTML = PALETTE.map(p =>
    `<button type="button" class="kit-ball" data-z="${p.Z}" style="background:${cpkColor(p.Z)}" title="${ELEMENTS[p.Z - 1].name}">${p.name}</button>`
  ).join('');

  const stage = root.querySelector('#kit-stage');
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  stage.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#0c1016');
  const camera = new THREE.PerspectiveCamera(40, 1, 0.05, 80);
  camera.position.set(4, -7, 5);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xffffff, 0x222222, 1.1));
  const key = new THREE.DirectionalLight(0xffffff, 0.7);
  key.position.set(3, -2, 6);
  scene.add(key);

  const plane = new THREE.Mesh(
    new THREE.PlaneGeometry(14, 10),
    new THREE.MeshStandardMaterial({ color: '#1a2030', roughness: 0.9 })
  );
  plane.rotation.x = 0;
  scene.add(plane);
  const grid = new THREE.GridHelper(14, 14, 0x334, 0x223);
  grid.rotation.x = Math.PI / 2;
  scene.add(grid);

  const mol = new THREE.Group();
  scene.add(mol);
  const cloud = new THREE.Group();
  scene.add(cloud);

  const state = { atoms: [], bonds: [], pick: 6, selected: -1, mode: 'place' };
  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  let dragging = -1;

  function setStatus(t) { root.querySelector('#kit-status').textContent = t; }

  function rebuild() {
    for (const c of [...mol.children]) {
      mol.remove(c);
      c.geometry?.dispose();
      c.material?.dispose();
    }
    state.atoms.forEach((a, i) => {
      const m = new THREE.Mesh(
        new THREE.SphereGeometry(0.28 + 0.0012 * covalentRadius(a.Z), 28, 20),
        new THREE.MeshStandardMaterial({
          color: cpkColor(a.Z), roughness: 0.28, metalness: 0.08,
          emissive: i === state.selected ? 0x224466 : 0,
        })
      );
      m.position.set(a.x, a.y, a.z);
      m.userData.i = i;
      mol.add(m);
    });
    state.bonds.forEach(b => {
      const A = state.atoms[b.a], B = state.atoms[b.b];
      if (!A || !B) return;
      const dir = new THREE.Vector3(B.x - A.x, B.y - A.y, B.z - A.z);
      const len = dir.length();
      const n = b.order;
      for (let k = 0; k < n; k++) {
        const stick = new THREE.Mesh(
          new THREE.CylinderGeometry(0.045, 0.045, len, 10),
          new THREE.MeshStandardMaterial({ color: '#e8eef6', roughness: 0.4 })
        );
        stick.position.set((A.x + B.x) / 2, (A.y + B.y) / 2, (A.z + B.z) / 2);
        stick.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
        const side = new THREE.Vector3().crossVectors(dir, new THREE.Vector3(0, 0, 1));
        if (side.length() < 1e-4) side.set(1, 0, 0);
        side.normalize().multiplyScalar((k - (n - 1) / 2) * 0.12);
        stick.position.add(side);
        stick.userData.bond = state.bonds.indexOf(b);
        mol.add(stick);
      }
    });
  }

  function autoBonds() {
    const kept = [];
    for (let i = 0; i < state.atoms.length; i++) {
      for (let j = 0; j < i; j++) {
        const A = state.atoms[i], B = state.atoms[j];
        const d = Math.hypot(A.x - B.x, A.y - B.y, A.z - B.z);
        const existing = state.bonds.find(b => (b.a === i && b.b === j) || (b.a === j && b.b === i));
        if (d < snapDistance(A.Z, B.Z)) kept.push(existing ?? { a: j, b: i, order: 1 });
      }
    }
    state.bonds = kept;
  }

  function place(x, y) {
    state.atoms.push({ Z: state.pick, x, y, z: 0 });
    autoBonds();
    rebuild();
    clearQuantum();
  }

  function clearQuantum() {
    for (const c of [...cloud.children]) { cloud.remove(c); c.geometry?.dispose(); c.material?.dispose(); }
    root.querySelector('#kit-info').innerHTML = '';
  }

  atomsBox.querySelectorAll('button').forEach(b => {
    b.addEventListener('click', () => {
      state.pick = Number(b.dataset.z);
      state.mode = 'place';
      atomsBox.querySelectorAll('button').forEach(x => x.classList.toggle('on', x === b));
    });
  });
  atomsBox.querySelector('[data-z="6"]').classList.add('on');

  root.querySelector('#kit-clear').onclick = () => { state.atoms = []; state.bonds = []; state.selected = -1; rebuild(); clearQuantum(); };
  root.querySelector('#kit-del').onclick = () => {
    if (state.selected < 0) return;
    const i = state.selected;
    state.atoms.splice(i, 1);
    state.bonds = state.bonds.filter(b => b.a !== i && b.b !== i).map(b => ({
      ...b, a: b.a > i ? b.a - 1 : b.a, b: b.b > i ? b.b - 1 : b.b,
    }));
    state.selected = -1;
    rebuild();
    clearQuantum();
  };
  root.querySelector('#kit-bond').onclick = () => {
    if (!state.bonds.length) return;
    const b = state.bonds[state.bonds.length - 1];
    b.order = b.order % 3 + 1;
    rebuild();
  };

  function eventPos(e) {
    const r = renderer.domElement.getBoundingClientRect();
    pointer.x = ((e.clientX - r.left) / r.width) * 2 - 1;
    pointer.y = -((e.clientY - r.top) / r.height) * 2 + 1;
  }

  renderer.domElement.addEventListener('pointerdown', (e) => {
    eventPos(e);
    raycaster.setFromCamera(pointer, camera);
    const hits = raycaster.intersectObjects(mol.children, false);
    const atomHit = hits.find(h => h.object.userData.i !== undefined);
    if (atomHit) {
      dragging = atomHit.object.userData.i;
      state.selected = dragging;
      controls.enabled = false;
      rebuild();
      return;
    }
    const planeHit = raycaster.intersectObject(plane)[0];
    if (planeHit && state.mode === 'place') place(planeHit.point.x, planeHit.point.y);
  });
  renderer.domElement.addEventListener('pointermove', (e) => {
    if (dragging < 0) return;
    eventPos(e);
    raycaster.setFromCamera(pointer, camera);
    const hit = raycaster.intersectObject(plane)[0];
    if (!hit) return;
    state.atoms[dragging].x = hit.point.x;
    state.atoms[dragging].y = hit.point.y;
    autoBonds();
    rebuild();
  });
  window.addEventListener('pointerup', () => { dragging = -1; controls.enabled = true; });

  async function calculate() {
    if (state.atoms.length < 1) { setStatus('Metti almeno un atomo.'); return; }
    if (state.atoms.length > 16) { setStatus('Massimo 16 atomi per il calcolo Hartree–Fock in questa vista.'); return; }
    setStatus('Ottimizzazione Hartree–Fock 6-31G…');
    const cx0 = state.atoms.reduce((s, a) => s + a.x, 0) / state.atoms.length;
    const cy0 = state.atoms.reduce((s, a) => s + a.y, 0) / state.atoms.length;
    const atoms = state.atoms.map(a => ({ Z: a.Z, xyz: [(a.x - cx0) * BOHR, (a.y - cy0) * BOHR, a.z * BOHR] }));
    try {
      const r = await request('optimize', { atoms, opts: { basis: '6-31G', maxSteps: 40 } });
      const opt = r.atoms;
      state.atoms = opt.map(a => ({ Z: a.Z, x: a.xyz[0] / BOHR, y: a.xyz[1] / BOHR, z: a.xyz[2] / BOHR }));
      autoBonds();
      rebuild();
      const analysis = analyzeStructure({ atoms: state.atoms.map(a => ({ Z: a.Z })), bonds: state.bonds, charge: 0 });
      const rows = analysis.atoms.map(a =>
        `<tr><td>${a.symbol}${a.index + 1}</td><td>${a.hybridization}</td><td>${a.vsepr?.molecularGeometry ?? '—'}</td><td>${a.formalCharge}</td><td>${a.oxidation}</td></tr>`
      ).join('');
      root.querySelector('#kit-info').innerHTML = `
        <p>E = ${r.summary.energy.toFixed(4)} Ha · gruppo ${r.summary.pointGroup} · gap ${r.summary.gapEV?.toFixed(2) ?? '—'} eV</p>
        <table class="kit-table"><thead><tr><th>Atomo</th><th>Ibridazione</th><th>VSEPR</th><th>Carica formale</th><th>n.o.</th></tr></thead><tbody>${rows}</tbody></table>
        <p class="kit-note">Ibridazione da domini VSEPR (legami + coppie). Gli orbitali sotto sono il risultato del calcolo, non un disegno.</p>`;
      setStatus('Densità e HOMO…');
      await showQuantum(opt, r.summary);
      setStatus('Fatto. I punti sono campioni di Born della densità HF.');
    } catch (err) {
      setStatus(`Calcolo non riuscito: ${err.message}`);
    }
  }

  async function showQuantum(optAtoms, summary) {
    clearQuantum();
    const center = [0, 0, 0];
    const half = 7;
    const dens = await request('grid', { kind: 'density', res: 28, half, center });
    const pts = sampleDensity(dens, 1800);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pts, 3));
    const mat = new THREE.PointsMaterial({ color: 0x7ec8ff, size: 0.06, transparent: true, opacity: 0.55, depthWrite: false });
    const points = new THREE.Points(geo, mat);
    cloud.add(points);
    const homo = summary.homo;
    if (homo !== undefined && homo >= 0) {
      const mo = await request('grid', { kind: 'mo', index: homo, spin: 'α', res: 28, half, center });
      const iso = 0.05;
      for (const sign of [1, -1]) {
        const mc = marchingCubes(mo, iso, sign);
        if (!mc?.positions?.length) continue;
        const g = new THREE.BufferGeometry();
        g.setAttribute('position', new THREE.Float32BufferAttribute(mc.positions, 3));
        if (mc.normals) g.setAttribute('normal', new THREE.Float32BufferAttribute(mc.normals, 3));
        const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({
          color: sign > 0 ? 0x3ecbff : 0xe07a5f, transparent: true, opacity: 0.35,
          side: THREE.DoubleSide, roughness: 0.5,
        }));
        m.scale.setScalar(1 / BOHR);
        cloud.add(m);
      }
    }
  }

  function sampleDensity(grid, n) {
    const { values, res, half, step } = grid;
    let max = 0;
    for (let i = 0; i < values.length; i++) if (values[i] > max) max = values[i];
    const out = [];
    let guard = 0;
    while (out.length < n * 3 && guard < n * 40) {
      guard++;
      const i = (Math.random() * values.length) | 0;
      if (Math.random() > values[i] / max) continue;
      const x = i % res, y = ((i / res) | 0) % res, z = (i / (res * res)) | 0;
      out.push((-half + x * step) / BOHR, (-half + y * step) / BOHR, (-half + z * step) / BOHR);
    }
    return out;
  }

  root.querySelector('#kit-calc').onclick = calculate;

  let raf = 0;
  const loop = () => {
    raf = requestAnimationFrame(loop);
    controls.update();
    renderer.render(scene, camera);
  };
  const resize = () => {
    const w = stage.clientWidth || 640, h = stage.clientHeight || 480;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  new ResizeObserver(resize).observe(stage);
  resize();
  loop();

  return {
    dispose() {
      cancelAnimationFrame(raf);
      cancelAll();
      renderer.dispose();
    },
  };
}
