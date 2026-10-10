// Scena 3D: nuvole di punti (densità di probabilità) e isosuperfici, assi in picometri.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { marchingCubes } from './marching.js';

const BOHR_PM = 52.917721090;

function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function dotTexture() {
  const c = document.createElement('canvas');
  c.width = c.height = 64;
  const g = c.getContext('2d');
  const grad = g.createRadialGradient(32, 32, 0, 32, 32, 32);
  grad.addColorStop(0, 'rgba(255,255,255,1)');
  grad.addColorStop(0.35, 'rgba(255,255,255,0.8)');
  grad.addColorStop(1, 'rgba(255,255,255,0)');
  g.fillStyle = grad;
  g.fillRect(0, 0, 64, 64);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

function niceStep(x) {
  const p = Math.pow(10, Math.floor(Math.log10(x)));
  const f = x / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
}

export function textSprite(text, color, sizePx = 28) {
  const c = document.createElement('canvas');
  const g = c.getContext('2d');
  const font = `500 ${sizePx}px "IBM Plex Mono", ui-monospace, monospace`;
  g.font = font;
  const w = Math.ceil(g.measureText(text).width) + 8;
  c.width = w;
  c.height = sizePx + 8;
  g.font = font;
  g.fillStyle = color;
  g.textBaseline = 'middle';
  g.fillText(text, 4, c.height / 2);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  const mat = new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true });
  const s = new THREE.Sprite(mat);
  s.userData.aspect = w / c.height;
  s.renderOrder = 10;
  return s;
}

export class Viewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.localClippingEnabled = true;
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(40, 1, 0.001, 1000);
    this.camera.up.set(0, 0, 1); // asse z verso l'alto, come nei testi di chimica
    this.camera.position.set(6, -9, 5);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.autoRotateSpeed = 1.2;

    this.scene.add(new THREE.HemisphereLight(0xf0f0f0, 0x606060, 0.9));
    const key = new THREE.DirectionalLight(0xffffff, 0.85);
    key.position.set(3, -4, 6);
    this.camera.add(key);
    this.scene.add(this.camera);

    this.content = new THREE.Group();
    this.scene.add(this.content);
    this.axes = new THREE.Group();
    this.scene.add(this.axes);
    this.nuclei = new THREE.Group();
    this.scene.add(this.nuclei);

    this.dot = dotTexture();
    this.clipPlane = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
    this.clipping = false;
    this.extent = 5;

    this.resize();
    new ResizeObserver(() => this.resize()).observe(container);
    this.frameCallbacks = new Set();
    this.clock = new THREE.Clock();
    this.renderer.setAnimationLoop(() => {
      const t = this.clock.getElapsedTime();
      for (const cb of this.frameCallbacks) cb(t);
      this.controls.update();
      this.renderer.render(this.scene, this.camera);
    });
    this.applyTheme();
  }

  resize() {
    const w = this.container.clientWidth;
    const h = this.container.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  applyTheme() {
    this.dark = cssVar('--scene-mode') !== 'light';
    this.scene.background = new THREE.Color(cssVar('--scene-bg') || '#080b10');
    this.content.traverse(o => {
      if (o.isPoints) {
        o.material.blending = this.dark ? THREE.AdditiveBlending : THREE.NormalBlending;
        o.material.needsUpdate = true;
      }
    });
    this.buildAxes(this.extent);
  }

  clear() {
    for (const child of [...this.content.children]) {
      this.content.remove(child);
      child.traverse?.(o => { o.geometry?.dispose(); o.material?.map?.dispose(); o.material?.dispose(); });
    }
    for (const child of [...this.nuclei.children]) this.nuclei.remove(child);
    this.frameCallbacks.clear();
  }

  /** Aggiunge un oggetto three.js qualsiasi al contenuto della scena. */
  add(object) {
    this.content.add(object);
    return object;
  }

  /** Mostra o nasconde gli assi cartesiani. */
  setAxesVisible(on) {
    this.axes.visible = on;
  }

  /** Superficie con colori per vertice (per esempio il potenziale elettrostatico). */
  addColoredSurface(positions, normals, colors, { opacity = 0.9 } = {}) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const mat = new THREE.MeshStandardMaterial({
      vertexColors: true, roughness: 0.45, metalness: 0.0,
      transparent: opacity < 1, opacity, side: THREE.DoubleSide, depthWrite: opacity >= 1,
      clippingPlanes: this.clipping ? [this.clipPlane] : [],
    });
    const mesh = new THREE.Mesh(geo, mat);
    this.content.add(mesh);
    return mesh;
  }

  setAutoRotate(on) {
    this.controls.autoRotate = on;
  }

  setClipping(on) {
    this.clipping = on;
    this.content.traverse(o => {
      if (o.material) {
        o.material.clippingPlanes = on ? [this.clipPlane] : [];
        o.material.needsUpdate = true;
      }
    });
  }

  /** Inquadra un oggetto di raggio `extent` (in a₀). */
  frame(extent, keepView = false) {
    const changed = Math.abs(extent - this.extent) / this.extent > 0.02;
    this.extent = extent;
    const dist = extent * 3.1;
    if (!keepView || changed) {
      const dir = this.camera.position.clone().sub(this.controls.target);
      if (dir.lengthSq() === 0) dir.set(0.55, -0.8, 0.45);
      dir.normalize().multiplyScalar(dist);
      this.controls.target.set(0, 0, 0);
      this.camera.position.copy(dir);
    }
    this.camera.near = extent / 500;
    this.camera.far = extent * 60;
    this.camera.updateProjectionMatrix();
    this.controls.minDistance = extent * 0.2;
    this.controls.maxDistance = extent * 20;
    this.buildAxes(extent);
  }

  resetView() {
    this.camera.position.set(0.55, -0.8, 0.45).normalize().multiplyScalar(this.extent * 3.1);
    this.controls.target.set(0, 0, 0);
  }

  /** Assi cartesiani con tacche a intervalli "rotondi" in picometri. */
  buildAxes(extent) {
    for (const c of [...this.axes.children]) {
      this.axes.remove(c);
      c.geometry?.dispose();
      c.material?.map?.dispose();
      c.material?.dispose();
    }
    const color = cssVar('--scene-axis') || '#5b6679';
    const textColor = cssVar('--scene-label') || '#9aa4b5';
    const L = extent * 1.15;
    const extentPm = L * BOHR_PM;
    const stepPm = niceStep(extentPm / 3);
    const step = stepPm / BOHR_PM;
    const verts = [];
    const tick = extent * 0.025;
    const dirs = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    dirs.forEach((d) => {
      verts.push(-L * d[0], -L * d[1], -L * d[2], L * d[0], L * d[1], L * d[2]);
      for (let t = step; t < L; t += step) {
        for (const s of [-1, 1]) {
          const p = [s * t * d[0], s * t * d[1], s * t * d[2]];
          // tacca perpendicolare all'asse
          const q = d[2] === 1 ? [tick, 0, 0] : [0, 0, tick];
          verts.push(p[0] - q[0], p[1] - q[1], p[2] - q[2], p[0] + q[0], p[1] + q[1], p[2] + q[2]);
        }
      }
    });
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
    const lines = new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.7 }));
    this.axes.add(lines);
    const labelH = extent * 0.09;
    ['x', 'y', 'z'].forEach((name, i) => {
      const s = textSprite(name, textColor, 36);
      const d = dirs[i];
      s.position.set(d[0] * L * 1.07, d[1] * L * 1.07, d[2] * L * 1.07);
      s.scale.set(labelH * s.userData.aspect, labelH, 1);
      this.axes.add(s);
    });
    const unit = textSprite(`${formatPm(stepPm)}`, textColor, 26);
    const uh = extent * 0.06;
    unit.position.set(step, 0, -extent * 0.07);
    unit.scale.set(uh * unit.userData.aspect, uh, 1);
    this.axes.add(unit);
    this.tickPm = stepPm;
  }

  addNucleus(position = [0, 0, 0], label = null) {
    const r = this.extent * 0.006;
    const mesh = new THREE.Mesh(
      new THREE.SphereGeometry(r, 20, 14),
      new THREE.MeshStandardMaterial({ color: cssVar('--nucleus') || '#c9a227', emissive: 0x000000, emissiveIntensity: 0, roughness: 0.6 }),
    );
    mesh.position.set(...position);
    this.nuclei.add(mesh);
    if (label) {
      const s = textSprite(label, cssVar('--scene-label') || '#9aa4b5', 30);
      const h = this.extent * 0.07;
      s.position.set(position[0] + r * 2.5, position[1], position[2] + r * 2.5);
      s.scale.set(h * s.userData.aspect, h, 1);
      this.nuclei.add(s);
    }
  }

  /** Nuvola di punti; `colors` Float32Array RGB per punto. */
  addPoints(positions, colors, { size = 1, opacity = 0.9 } = {}) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const mat = new THREE.PointsMaterial({
      size: this.extent * 0.016 * size,
      map: this.dot,
      vertexColors: true,
      transparent: true,
      opacity,
      depthWrite: false,
      blending: THREE.NormalBlending,
      clippingPlanes: this.clipping ? [this.clipPlane] : [],
    });
    const pts = new THREE.Points(geo, mat);
    this.content.add(pts);
    return pts;
  }

  /** Isosuperficie dal campo campionato su griglia. */
  addSurface(gridField, iso, sign, color, { opacity = 0.82 } = {}) {
    const { positions, normals } = marchingCubes(gridField, iso, sign);
    if (positions.length === 0) return null;
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
    const mat = new THREE.MeshStandardMaterial({
      color,
      roughness: 0.38,
      metalness: 0.05,
      transparent: opacity < 1,
      opacity,
      side: THREE.DoubleSide,
      depthWrite: opacity >= 1,
      clippingPlanes: this.clipping ? [this.clipPlane] : [],
    });
    const mesh = new THREE.Mesh(geo, mat);
    this.content.add(mesh);
    return mesh;
  }

  snapshot() {
    return this.renderer.domElement.toDataURL('image/png');
  }
}

export function formatPm(pm) {
  if (pm >= 1) return `${+pm.toPrecision(3)} pm`;
  return `${+(pm * 1000).toPrecision(3)} fm`;
}

export { BOHR_PM };
