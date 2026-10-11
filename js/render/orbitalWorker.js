import { prepareOrbitalScene } from './orbitalScene.js';

let latestId = 0;
self.onmessage = async ({ data: { id, payload } }) => {
  latestId = id;
  try {
    const scene = await prepareOrbitalScene(payload, () => id !== latestId);
    const transfer = scene.surfaces.flatMap(s => [s.positions.buffer, s.normals.buffer]);
    if (scene.cloud) transfer.push(scene.cloud.positions.buffer, scene.cloud.groups.buffer);
    self.postMessage({ id, scene }, transfer);
  } catch (error) {
    self.postMessage({ id, error: error.message });
  }
};
