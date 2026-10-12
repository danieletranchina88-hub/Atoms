// Historical models only; the live quantum simulation uses its unchanged renderer.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { rutherford } from '../course/chemistry.js';

export function makeStage(container, { getPaused = () => false } = {}) {
  const wrap = document.createElement('div');
  wrap.className = 'course-model-stage';
  container.replaceChildren(wrap);
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  wrap.appendChild(renderer.domElement);
  renderer.domElement.setAttribute(
    'aria-label',
    'Modello atomico storico 3D, trascina per ruotare',
  );
  const scene = new THREE.Scene(),
    camera = new THREE.PerspectiveCamera(40, 1, 0.01, 200);
  scene.background = new THREE.Color('#101821');
  camera.position.set(0, -8, 4);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  scene.add(new THREE.HemisphereLight(0xffffff, 0x304050, 2));
  const light = new THREE.DirectionalLight(0xffffff, 2);
  light.position.set(3, -4, 5);
  scene.add(light);
  const content = new THREE.Group();
  scene.add(content);
  let disposed = false,
    previous = 0,
    time = 0,
    paused = false;
  const reduced = matchMedia('(prefers-reduced-motion: reduce)');
  const resize = () => {
    const w = wrap.clientWidth || 600,
      h = wrap.clientHeight || 320;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  };
  const ro = new ResizeObserver(resize);
  ro.observe(wrap);
  resize();
  const clear = () => {
    content.userData.tick = null;
    content.rotation.set(0, 0, 0);
    time = 0;
    const materials = new Set(),
      geometries = new Set();
    for (const child of [...content.children]) {
      child.traverse((o) => {
        if (o.geometry) geometries.add(o.geometry);
        if (o.material)
          (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) =>
            materials.add(m),
          );
      });
      content.remove(child);
    }
    geometries.forEach((g) => g.dispose());
    materials.forEach((m) => m.dispose());
  };
  renderer.setAnimationLoop((now) => {
    if (disposed) return;
    const dt = previous ? Math.min((now - previous) / 1000, 0.05) : 0;
    previous = now;
    if (document.hidden || !wrap.isConnected) return;
    if (!paused && !getPaused() && !reduced.matches) time += dt;
    content.userData.tick?.(time);
    controls.update();
    renderer.render(scene, camera);
  });
  return {
    scene,
    camera,
    content,
    controls,
    clear,
    setPaused: (on) => {
      paused = on;
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      renderer.setAnimationLoop(null);
      ro.disconnect();
      controls.dispose();
      clear();
      renderer.dispose();
      renderer.forceContextLoss();
      wrap.remove();
    },
  };
}
function sphere(radius, color, opacity = 1) {
  return new THREE.Mesh(
    new THREE.SphereGeometry(radius, 32, 24),
    new THREE.MeshStandardMaterial({
      color,
      roughness: 0.36,
      transparent: opacity < 1,
      opacity,
      depthWrite: opacity === 1,
    }),
  );
}
export function showDalton(stage) {
  stage.clear();
  const a = sphere(0.7, '#f0b429'),
    b = sphere(0.55, '#5ab8e7');
  a.position.x = -1;
  b.position.x = 1;
  stage.content.add(a, b);
  stage.controls.target.set(0, 0, 0);
  stage.camera.position.set(0, -8, 3);
}
export function showThomson(stage) {
  stage.clear();
  stage.content.add(sphere(1.4, '#f0b429', 0.22));
  // Static particles in the historical positive background: no fabricated oscillation.
  for (let i = 0; i < 8; i++) {
    const a = sphere(0.11, '#5ab8e7');
    const z = 1 - (2 * (i + 0.5)) / 8,
      r = Math.sqrt(1 - z * z),
      phi = i * 2.39996323;
    a.position.set(
      0.95 * r * Math.cos(phi),
      0.95 * r * Math.sin(phi),
      0.95 * z,
    );
    stage.content.add(a);
  }
  stage.controls.target.set(0, 0, 0);
  stage.camera.position.set(0, -7, 3);
}
export function showRutherford(stage, params = {}) {
  stage.clear();
  const model = rutherford(params),
    radius = Math.max(...model.path.map((p) => Math.hypot(p.x, p.y))),
    scale = 4 / radius;
  const points = model.path.map(
    (p) => new THREE.Vector3(p.x * scale, p.y * scale, 0),
  );
  stage.content.add(
    sphere(0.07, '#f58b76'),
    new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(points),
      new THREE.LineBasicMaterial({ color: '#5ab8e7' }),
    ),
  );
  const alpha = sphere(0.05, '#f0b429');
  stage.content.add(alpha);
  stage.content.userData.tick = (time) => {
    const u = (time / 6) % 1;
    let i = model.path.findIndex((p) => p.t >= u);
    if (i < 1) i = 1;
    const f =
      (u - model.path[i - 1].t) / (model.path[i].t - model.path[i - 1].t);
    alpha.position.lerpVectors(points[i - 1], points[i], f);
  };
  stage.controls.target.set(0, 0, 0);
  stage.camera.position.set(0, -0.1, 11);
  return model;
}

export function showIsomers(stage) {
  stage.clear();
  const build = (atoms, bonds, offset) => {
    const group = new THREE.Group();
    group.position.x = offset;
    atoms.forEach(([el, x, y, z]) => {
      const atom = sphere(
        el === 'H' ? 0.12 : 0.22,
        el === 'O' ? '#f58b76' : el === 'H' ? '#e8edf2' : '#8393a3',
      );
      atom.position.set(x, y, z);
      atom.userData.element = el;
      group.add(atom);
    });
    bonds.forEach(([i, j]) => {
      const a = new THREE.Vector3(...atoms[i].slice(1)),
        b = new THREE.Vector3(...atoms[j].slice(1)),
        delta = b.clone().sub(a);
      const mesh = new THREE.Mesh(
        new THREE.CylinderGeometry(0.035, 0.035, delta.length(), 8),
        new THREE.MeshStandardMaterial({ color: '#647789', roughness: 0.5 }),
      );
      mesh.position.copy(a).add(b).multiplyScalar(0.5);
      mesh.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        delta.normalize(),
      );
      group.add(mesh);
    });
    stage.content.add(group);
  };
  build(
    [
      ['C', -1.4, 0, 0],
      ['C', 0.1, 0, 0],
      ['O', 0.85, 1.2, 0],
      ['H', -1.8, 0.95, 0],
      ['H', -1.8, -0.45, 0.85],
      ['H', -1.8, -0.45, -0.85],
      ['H', 0.55, -0.45, 0.85],
      ['H', 0.55, -0.45, -0.85],
      ['H', 1.7, 1.2, 0],
    ],
    [
      [0, 1],
      [1, 2],
      [0, 3],
      [0, 4],
      [0, 5],
      [1, 6],
      [1, 7],
      [2, 8],
    ],
    -2.3,
  );
  build(
    [
      ['O', 0, 0.5, 0],
      ['C', -1.2, -0.35, 0],
      ['C', 1.2, -0.35, 0],
      ['H', -1.7, 0.45, 0],
      ['H', -1.55, -0.9, 0.85],
      ['H', -1.55, -0.9, -0.85],
      ['H', 1.7, 0.45, 0],
      ['H', 1.55, -0.9, 0.85],
      ['H', 1.55, -0.9, -0.85],
    ],
    [
      [0, 1],
      [0, 2],
      [1, 3],
      [1, 4],
      [1, 5],
      [2, 6],
      [2, 7],
      [2, 8],
    ],
    2.3,
  );
  stage.controls.target.set(0, 0, 0);
  stage.camera.position.set(0, -13, 5);
}
