// Scene 3D del corso. Ogni modello è etichettato con chi lo ha proposto.
// Dalton: sfere indivisibili. Thomson: panettone. Rutherford: nucleo e scattering coulombiano.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

export function makeStage(container) {
  container.innerHTML = '';
  const wrap = document.createElement('div');
  wrap.className = 'course-stage';
  container.appendChild(wrap);
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  wrap.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#07090d');
  const camera = new THREE.PerspectiveCamera(40, 1, 0.05, 200);
  camera.position.set(6, -8, 4);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xf4f4f4, 0x303030, 1));
  const key = new THREE.DirectionalLight(0xffffff, 0.8);
  key.position.set(4, -3, 6);
  scene.add(key);
  const content = new THREE.Group();
  scene.add(content);
  const clock = new THREE.Clock();
  let raf = 0;
  const resize = () => {
    const w = wrap.clientWidth || 640, h = 280;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  new ResizeObserver(resize).observe(wrap);
  resize();
  const loop = () => {
    raf = requestAnimationFrame(loop);
    const t = clock.getElapsedTime();
    content.userData.tick?.(t);
    controls.update();
    renderer.render(scene, camera);
  };
  loop();
  return {
    scene, camera, content, controls,
    clear() {
      for (const c of [...content.children]) {
        content.remove(c);
        c.traverse?.(o => { o.geometry?.dispose(); o.material?.dispose(); });
      }
      content.userData.tick = null;
    },
    dispose() { cancelAnimationFrame(raf); renderer.dispose(); wrap.remove(); },
  };
}

function sphere(r, color, emissive = 0) {
  const m = new THREE.MeshStandardMaterial({ color, roughness: 0.45, metalness: 0, emissive, emissiveIntensity: emissive ? 0.4 : 0 });
  return new THREE.Mesh(new THREE.SphereGeometry(r, 32, 24), m);
}

export function showDalton(stage) {
  stage.clear();
  const a = sphere(0.7, '#c9a227');
  const b = sphere(0.55, '#e07a5f');
  b.position.set(1.6, 0, 0);
  stage.content.add(a, b);
  stage.content.userData.tick = (t) => { stage.content.rotation.z = t * 0.15; };
  stage.controls.target.set(0.8, 0, 0);
}

export function showThomson(stage) {
  stage.clear();
  const pudding = new THREE.Mesh(
    new THREE.SphereGeometry(1.4, 48, 36),
    new THREE.MeshStandardMaterial({ color: '#f2c14e', transparent: true, opacity: 0.35, roughness: 0.6 })
  );
  stage.content.add(pudding);
  const electrons = [];
  for (let i = 0; i < 8; i++) {
    const e = sphere(0.12, '#3ecbff', 0x2266aa);
    const phi = Math.acos(2 * Math.random() - 1);
    const th = Math.random() * Math.PI * 2;
    const r = 0.4 + Math.random() * 0.8;
    e.position.set(r * Math.sin(phi) * Math.cos(th), r * Math.sin(phi) * Math.sin(th), r * Math.cos(phi));
    e.userData.base = e.position.clone();
    electrons.push(e);
    stage.content.add(e);
  }
  stage.content.userData.tick = (t) => {
    electrons.forEach((e, i) => {
      e.position.copy(e.userData.base);
      e.position.x += Math.sin(t * 1.4 + i) * 0.06;
      e.position.y += Math.cos(t * 1.1 + i) * 0.06;
    });
  };
}

export function showRutherford(stage, { Z = 79, E = 5, b = 20 } = {}) {
  stage.clear();
  const nucleus = sphere(0.18, '#e07a5f', 0x882222);
  stage.content.add(nucleus);
  const k = 1.43996;
  const theta = 2 * Math.atan((Z * 2 * k) / (2 * E * b));
  const curve = new THREE.CatmullRomCurve3(
    Array.from({ length: 40 }, (_, i) => {
      const s = (i / 39) * 2 - 1;
      const along = s * 4;
      const off = (b / 40) * (1 - Math.abs(s)) + (s > 0 ? Math.tan(theta / 2) * along * 0.25 : 0);
      return new THREE.Vector3(along, off, 0);
    })
  );
  const geo = new THREE.BufferGeometry().setFromPoints(curve.getPoints(80));
  stage.content.add(new THREE.Line(geo, new THREE.LineBasicMaterial({ color: '#3ecbff' })));
  const alpha = sphere(0.1, '#c9a227', 0x886622);
  stage.content.add(alpha);
  stage.content.userData.tick = (t) => {
    const u = (t * 0.25) % 1;
    alpha.position.copy(curve.getPoint(u));
  };
  stage.controls.target.set(0, 0.4, 0);
}

const CPK = { H: '#ffffff', C: '#909090', O: '#ff0d0d', N: '#3050f8', Fe: '#e06633', Al: '#bfa6a6', S: '#ffff30', Cl: '#1ff01f' };

function ball(el, r = 0.28) {
  return sphere(el === 'H' ? 0.18 : r, CPK[el] || '#c9a227');
}

function molecule(formula) {
  const g = new THREE.Group();
  const f = formula.replace(/\s/g, '');
  if (f === 'H2O') {
    const o = ball('O'); const h1 = ball('H'); const h2 = ball('H');
    h1.position.set(0.7, 0.4, 0); h2.position.set(-0.7, 0.4, 0);
    g.add(o, h1, h2);
  } else if (f === 'CH4') {
    const c = ball('C');
    [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].forEach(p => {
      const h = ball('H'); h.position.set(p[0] * 0.55, p[1] * 0.55, p[2] * 0.55); g.add(h);
    });
    g.add(c);
  } else if (f === 'O2') {
    const a = ball('O'), b = ball('O'); b.position.x = 0.7; g.add(a, b);
  } else if (f === 'CO2') {
    const c = ball('C'), o1 = ball('O'), o2 = ball('O');
    o1.position.x = -0.8; o2.position.x = 0.8; g.add(c, o1, o2);
  } else if (f === 'Fe2O3' || f === 'Al2O3') {
    const metal = f.startsWith('Fe') ? 'Fe' : 'Al';
    const m1 = ball(metal), m2 = ball(metal);
    m1.position.set(-0.4, 0.3, 0); m2.position.set(0.4, 0.3, 0);
    [[0, -0.3, 0.4], [0, -0.3, -0.4], [0, 0.6, 0]].forEach(p => {
      const o = ball('O', 0.22); o.position.set(...p); g.add(o);
    });
    g.add(m1, m2);
  } else {
    g.add(ball(f.replace(/\d/g, '') || 'C'));
  }
  return g;
}

export function showReaction(stage, left, right, balanced) {
  stage.clear();
  let x = -3;
  for (const sp of left) {
    const n = Math.min(sp.coeff || 0, 4);
    for (let i = 0; i < n; i++) {
      const m = molecule(sp.formula);
      m.position.set(x, (i - n / 2) * 0.9, 0);
      stage.content.add(m);
    }
    x += 1.6;
  }
  x = 1.2;
  for (const sp of right) {
    const n = Math.min(sp.coeff || 0, 4);
    for (let i = 0; i < n; i++) {
      const m = molecule(sp.formula);
      m.position.set(x, (i - n / 2) * 0.9, 0);
      stage.content.add(m);
    }
    x += 1.6;
  }
  stage.content.userData.tick = (t) => { stage.content.rotation.z = Math.sin(t * 0.4) * 0.08; };
  stage.controls.target.set(0, 0, 0);
  stage.camera.position.set(0, -10, 4);
}

export function showElectronTransfer(stage, n) {
  stage.clear();
  const donor = sphere(0.55, '#e07a5f');
  const acceptor = sphere(0.55, '#3ecbff');
  donor.position.x = -2; acceptor.position.x = 2;
  stage.content.add(donor, acceptor);
  const electrons = [];
  const count = Math.min(n || 1, 8);
  for (let i = 0; i < count; i++) {
    const e = sphere(0.1, '#c9a227', 0x886622);
    e.userData.phase = i / count;
    electrons.push(e);
    stage.content.add(e);
  }
  stage.content.userData.tick = (t) => {
    electrons.forEach(e => {
      const u = (t * 0.35 + e.userData.phase) % 1;
      e.position.set(-2 + u * 4, Math.sin(u * Math.PI) * 0.8, 0);
    });
  };
}
