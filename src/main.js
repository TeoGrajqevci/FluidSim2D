// main.js
import { canvas, gl, simHeight, cScale, simWidth } from './common.js';
import { FlipFluid } from './flipFluid.js';

const scene = {
  gravity: -9.81,
  dt: 1.0 / 60.0,
  flipRatio: 0.7,
  numParticleIters: 2,
  numPressureIters: 70,
  overRelaxation: 1.9,
  compensateDrift: true,
  separateParticles: true,
  obstacleRadius: 0.1,
  obstacleVelX: 0.0,
  obstacleVelY: 0.0,
  showGrid: true,
  fluid: null,
  frameNr: 0,
  paused: false
};

function setupScene() {
  scene.obstacleRadius = 0.1;
  const res = 110;
  const tankHeight = simHeight, tankWidth = simWidth;
  const h = tankHeight / res, density = 800.0;
  const relWaterHeight = 0.7, relWaterWidth = 0.7;
  const r = 0.3 * h, dx = 2.0 * r, dy = Math.sqrt(3.0) / 2.0 * dx;
  const numX = Math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
  const numY = Math.floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
  const maxParticles = numX * numY;
  scene.fluid = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);
  scene.fluid.numParticles = numX * numY;
  let p = 0;
  for (let i = 0; i < numX; i++) {
    for (let j = 0; j < numY; j++) {
      scene.fluid.particlePos[p++] = h + r + dx * i + (j % 2 === 0 ? 0.0 : r);
      scene.fluid.particlePos[p++] = h + r + dy * j;
    }
  }
  const n = scene.fluid.fNumY;
  for (let i = 0; i < scene.fluid.fNumX; i++) {
    for (let j = 0; j < scene.fluid.fNumY; j++) {
      scene.fluid.s[i * n + j] = (i === 0 || i === scene.fluid.fNumX - 1 || j === 0) ? 0.0 : 1.0;
    }
  }
  setObstacle(1000.0, 0.0, true);
}

// Shaders and drawing
const pointVertexShader = `
  attribute vec2 attrPosition;
  attribute vec3 attrColor;
  uniform vec2 domainSize;
  uniform float pointSize;
  uniform float drawDisk;
  varying vec3 fragColor;
  varying float fragDrawDisk;
  void main() {
    vec4 transform = vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
    gl_Position = vec4(attrPosition * transform.xy + transform.zw, 0.0, 1.0);
    gl_PointSize = pointSize;
    fragColor = attrColor;
    fragDrawDisk = drawDisk;
  }
`;
const pointFragmentShader = `
  precision mediump float;
  varying vec3 fragColor;
  varying float fragDrawDisk;
  void main() {
    if (fragDrawDisk == 1.0) {
      float rx = 0.5 - gl_PointCoord.x;
      float ry = 0.5 - gl_PointCoord.y;
      if (rx * rx + ry * ry > 0.25) discard;
    }
    gl_FragColor = vec4(fragColor, 1.0);
  }
`;
function createShader(gl, vsSource, fsSource) {
  const vs = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vs, vsSource);
  gl.compileShader(vs);
  const fs = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fs, fsSource);
  gl.compileShader(fs);
  const shader = gl.createProgram();
  gl.attachShader(shader, vs);
  gl.attachShader(shader, fs);
  gl.linkProgram(shader);
  return shader;
}
let pointShader = null, gridVertBuffer = null, gridColorBuffer = null;
function draw() {
  gl.clearColor(0, 0, 0, 1);
  gl.clear(gl.COLOR_BUFFER_BIT);
  gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
  if (!pointShader)
    pointShader = createShader(gl, pointVertexShader, pointFragmentShader);
  if (!gridVertBuffer) {
    const f = scene.fluid;
    gridVertBuffer = gl.createBuffer();
    const cellCenters = new Float32Array(2 * f.fNumCells);
    let p = 0;
    for (let i = 0; i < f.fNumX; i++) {
      for (let j = 0; j < f.fNumY; j++) {
        cellCenters[p++] = (i + 0.5) * f.h;
        cellCenters[p++] = (j + 0.5) * f.h;
      }
    }
    gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, cellCenters, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }
  if (!gridColorBuffer)
    gridColorBuffer = gl.createBuffer();
  if (scene.showGrid) {
    const pointSize = 0.5 * scene.fluid.h / simWidth * canvas.width;
    gl.useProgram(pointShader);
    gl.uniform2f(gl.getUniformLocation(pointShader, 'domainSize'), simWidth, simHeight);
    gl.uniform1f(gl.getUniformLocation(pointShader, 'pointSize'), pointSize);
    gl.uniform1f(gl.getUniformLocation(pointShader, 'drawDisk'), 0.0);
    gl.bindBuffer(gl.ARRAY_BUFFER, gridVertBuffer);
    const posLoc = gl.getAttribLocation(pointShader, 'attrPosition');
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, gridColorBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, scene.fluid.cellColor, gl.DYNAMIC_DRAW);
    const colorLoc = gl.getAttribLocation(pointShader, 'attrColor');
    gl.enableVertexAttribArray(colorLoc);
    gl.vertexAttribPointer(colorLoc, 3, gl.FLOAT, false, 0, 0);
    gl.drawArrays(gl.POINTS, 0, scene.fluid.fNumCells);
    gl.disableVertexAttribArray(posLoc);
    gl.disableVertexAttribArray(colorLoc);
    gl.bindBuffer(gl.ARRAY_BUFFER, null);
  }
}

function setObstacle(x, y, reset) {
  let vx = 0.0, vy = 0.0;
  if (!reset) {
    vx = (x - scene.obstacleX) / scene.dt;
    vy = (y - scene.obstacleY) / scene.dt;
  }
  scene.obstacleX = x;
  scene.obstacleY = y;
  const r = scene.obstacleRadius;
  const f = scene.fluid, n = f.fNumY;
  for (let i = 1; i < f.fNumX - 2; i++) {
    for (let j = 1; j < f.fNumY - 2; j++) {
      f.s[i * n + j] = 1.0;
      let dx = (i + 0.5) * f.h - x,
          dy = (j + 0.5) * f.h - y;
      if (dx * dx + dy * dy < r * r) {
        f.s[i * n + j] = 0.0;
        f.u[i * n + j] = vx;
        f.u[(i + 1) * n + j] = vx;
        f.v[i * n + j] = vy;
        f.v[i * n + j + 1] = vy;
      }
    }
  }
  scene.obstacleVelX = vx;
  scene.obstacleVelY = vy;
}

let mouseDown = false;
function startDrag(x, y) {
  const bounds = canvas.getBoundingClientRect();
  let mx = x - bounds.left, my = y - bounds.top;
  mouseDown = true;
  x = mx / cScale;
  y = (canvas.height - my) / cScale;
  setObstacle(x, y, true);
  scene.paused = false;
}
function drag(x, y) {
  if (mouseDown) {
    const bounds = canvas.getBoundingClientRect();
    let mx = x - bounds.left, my = y - bounds.top;
    x = mx / cScale;
    y = (canvas.height - my) / cScale;
    setObstacle(x, y, false);
  }
}
function endDrag() {
  mouseDown = false;
  scene.obstacleVelX = 0.0;
  scene.obstacleVelY = 0.0;
}
canvas.addEventListener('mousedown', e => startDrag(e.x, e.y));
canvas.addEventListener('mouseup', endDrag);
canvas.addEventListener('mousemove', e => drag(e.x, e.y));
canvas.addEventListener('touchstart', e => startDrag(e.touches[0].clientX, e.touches[0].clientY));
canvas.addEventListener('touchend', endDrag);
canvas.addEventListener('touchmove', e => { e.preventDefault(); drag(e.touches[0].clientX, e.touches[0].clientY); }, {passive: false});

function simulate() {
  if (!scene.paused) {
    scene.fluid.simulate(
      scene.dt, scene.gravity, scene.flipRatio,
      scene.numPressureIters, scene.numParticleIters,
      scene.overRelaxation, scene.compensateDrift,
      scene.separateParticles,
      scene.obstacleX, scene.obstacleY, scene.obstacleRadius,
      scene.obstacleVelX, scene.obstacleVelY
    );
    scene.frameNr++;
  }
}

function update() {
  simulate();
  draw();
  requestAnimationFrame(update);
}

setupScene();
update();
