// Atomic-table-only rendering. The sandbox retains the original Viewer behavior.
import * as THREE from 'three';
import { Viewer } from './viewer.js';
import { marchingCubes } from './marching.js';

export class AtomViewer extends Viewer {
  constructor(container) {
    super(container);
    this.quantumViewActive = location.hash !== '#sandbox';
    this.quantumMaterials = new Set();
    this.quantumAnimated = !window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    this.quantumSpeed = 1;
    this.quantumTime = 0;
    let previousTime = 0;
    this.renderer.setAnimationLoop(() => {
      const t = this.clock.getElapsedTime();
      const dt = Math.min(Math.max(t - previousTime, 0), 0.05);
      previousTime = t;
      if (this.quantumViewActive) {
        if (document.hidden) return;
        if (this.quantumAnimated) this.quantumTime += dt * this.quantumSpeed;
        for (const mat of this.quantumMaterials) mat.uniforms.uTime.value = this.quantumTime;
      }
      // Same callback time, controls and drawing as Viewer when sandbox is active.
      for (const cb of this.frameCallbacks) cb(t);
      this.controls.update();
      this.renderer.render(this.scene, this.camera);
    });
  }

  setQuantumViewActive(on) { this.quantumViewActive = on; }

  setQuantumAnimation(on, speed = this.quantumSpeed) {
    this.quantumAnimated = on;
    this.quantumSpeed = speed;
  }

  clear() {
    if (!this.quantumViewActive) {
      super.clear();
      this.quantumMaterials.clear();
      this.quantumTime = 0;
      return;
    }
    for (const child of [...this.content.children]) {
      this.content.remove(child);
      child.traverse?.(o => { o.geometry?.dispose(); if (o.material?.map !== this.dot) o.material?.map?.dispose(); o.material?.dispose(); });
    }
    for (const child of [...this.nuclei.children]) {
      this.nuclei.remove(child);
      child.traverse(o => { o.geometry?.dispose(); o.material?.map?.dispose(); o.material?.dispose(); });
    }
    this.quantumMaterials.clear();
    this.quantumTime = 0;
    this.frameCallbacks.clear();
  }

  addGlowPoints(positions, colors, { size = 1, opacity = 0.9, animated = false, visibility = null } = {}) {
    if (!this.quantumViewActive) return super.addGlowPoints(positions, colors, { size, opacity });
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    const seeds = new Float32Array(positions.length / 3);
    for (let i = 0; i < seeds.length; i++) seeds[i] = ((i * 0.61803398875) % 1) * Math.PI * 2;
    geo.setAttribute('aSeed', new THREE.BufferAttribute(seeds, 1));
    geo.setAttribute('aVisibility', new THREE.BufferAttribute(visibility ?? new Float32Array(seeds.length), 1));
    const mat = new THREE.ShaderMaterial({
      transparent: true,
      depthWrite: false,
      blending: this.dark ? THREE.AdditiveBlending : THREE.NormalBlending,
      clipping: true,
      clippingPlanes: this.clipping ? [this.clipPlane] : [],
      uniforms: {
        uSize: { value: this.extent * 0.0085 * size },
        uOpacity: { value: opacity * Math.min(4, 60000 / Math.max(seeds.length, 1)) * 0.12 },
        uTime: { value: this.quantumTime },
        uScintillation: { value: animated && !visibility ? 1 : 0 },
        uCoherent: { value: visibility ? 1 : 0 },
        uViewport: { value: this.container.clientHeight * this.renderer.getPixelRatio() },
      },
      vertexShader: `
        #include <clipping_planes_pars_vertex>
        attribute vec3 color;
        attribute float aSeed;
        attribute float aVisibility;
        varying vec3 vColor;
        varying float vSeed;
        varying float vVisibility;
        uniform float uSize;
        uniform float uViewport;
        void main() {
          vColor = color;
          vSeed = aSeed;
          vVisibility = aVisibility;
          vec4 mv = modelViewMatrix * vec4(position, 1.0);
          gl_Position = projectionMatrix * mv;
          gl_PointSize = clamp(uSize * uViewport * projectionMatrix[1][1] / max(-mv.z, 0.0001), 1.5, 12.0);
          vec4 mvPosition = mv;
          #include <clipping_planes_vertex>
        }
      `,
      fragmentShader: `
        #include <clipping_planes_pars_fragment>
        varying vec3 vColor;
        varying float vSeed;
        varying float vVisibility;
        uniform float uOpacity;
        uniform float uTime;
        uniform float uScintillation;
        uniform float uCoherent;
        void main() {
          #include <clipping_planes_fragment>
          vec2 p = gl_PointCoord * 2.0 - 1.0;
          float r2 = dot(p, p);
          if (r2 > 1.0) discard;
          float core = exp(-r2 * 7.5);
          float halo = exp(-r2 * 2.0);
          // Random sample phases are independent of position: stationary density is unchanged.
          float scintillation = mix(1.0, 1.0 + 0.8 * sin(uTime * 2.4 + vSeed), uScintillation);
          // H 1s+2pz: ρ/q=1+2φ1φ2/(φ1²+φ2²) cos(ΔE t/ℏ).
          // One screen cycle is 4 seconds; the physical period is 405.3 attoseconds.
          float weight = 1.0 + uCoherent * vVisibility * cos(uTime * 1.57079632679);
          float a = uOpacity * (core * 0.8 + halo * 0.22) * scintillation * weight;
          gl_FragColor = vec4(vColor * (0.65 + core * 0.35), a);
          #include <tonemapping_fragment>
          #include <colorspace_fragment>
        }
      `,
    });
    const pts = new THREE.Points(geo, mat);
    if (animated || visibility) this.quantumMaterials.add(mat);
    pts.onBeforeRender = () => { mat.uniforms.uViewport.value = this.container.clientHeight * this.renderer.getPixelRatio(); };
    this.content.add(pts);
    return pts;
  }

  addSurface(gridField, iso, sign, color, options = {}) {
    if (!this.quantumViewActive) return super.addSurface(gridField, iso, sign, color, options);
    const { positions, normals } = marchingCubes(gridField, iso, sign);
    return this.addOrbitalSurface(positions, normals, color, options);
  }

  addOrbitalSurface(positions, normals, color, { opacity = 0.3, emissive = 0.15 } = {}) {
    if (positions.length === 0) return null;
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
    const mat = new THREE.MeshStandardMaterial({
      color,
      emissive: color,
      emissiveIntensity: emissive,
      roughness: 0.3,
      metalness: 0.04,
      transparent: opacity < 1,
      opacity,
      side: THREE.DoubleSide,
      depthWrite: opacity >= 1,
      clippingPlanes: this.clipping ? [this.clipPlane] : [],
    });
    mat.onBeforeCompile = shader => {
      shader.fragmentShader = shader.fragmentShader.replace('#include <opaque_fragment>', `
        float rim = pow(1.0 - abs(dot(normalize(normal), normalize(vViewPosition))), 2.0);
        outgoingLight += diffuse * rim * 0.35;
        diffuseColor.a *= 0.35 + 0.65 * rim;
        #include <opaque_fragment>
      `);
    };
    const mesh = new THREE.Mesh(geo, mat);
    this.content.add(mesh);
    return mesh;
  }

}
