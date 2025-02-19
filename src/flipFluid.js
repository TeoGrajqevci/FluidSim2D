// flipFluid.js
import { clamp, FLUID_CELL, AIR_CELL, SOLID_CELL } from './common.js';

export class FlipFluid {
  constructor(density, width, height, spacing, particleRadius, maxParticles) {
    this.density = density;
    this.fNumX = Math.floor(width / spacing) + 1;
    this.fNumY = Math.floor(height / spacing) + 1;
    this.h = Math.max(width / this.fNumX, height / this.fNumY);
    this.fInvSpacing = 1.0 / this.h;
    this.fNumCells = this.fNumX * this.fNumY;

    this.u = new Float32Array(this.fNumCells);
    this.v = new Float32Array(this.fNumCells);
    this.du = new Float32Array(this.fNumCells);
    this.dv = new Float32Array(this.fNumCells);
    this.prevU = new Float32Array(this.fNumCells);
    this.prevV = new Float32Array(this.fNumCells);
    this.p = new Float32Array(this.fNumCells);
    this.s = new Float32Array(this.fNumCells);
    this.cellType = new Int32Array(this.fNumCells);
    this.cellColor = new Float32Array(3 * this.fNumCells);

    this.maxParticles = maxParticles;
    this.particlePos = new Float32Array(2 * this.maxParticles);
    this.particleColor = new Float32Array(3 * this.maxParticles);
    for (let i = 0; i < this.maxParticles; i++) {
      this.particleColor[3 * i + 2] = 1.0;
    }
    this.particleVel = new Float32Array(2 * this.maxParticles);
    this.particleDensity = new Float32Array(this.fNumCells);
    this.particleRestDensity = 0.0;
    this.particleRadius = particleRadius;
    this.pInvSpacing = 1.0 / (2.2 * particleRadius);
    this.pNumX = Math.floor(width * this.pInvSpacing) + 1;
    this.pNumY = Math.floor(height * this.pInvSpacing) + 1;
    this.pNumCells = this.pNumX * this.pNumY;
    this.numCellParticles = new Int32Array(this.pNumCells);
    this.firstCellParticle = new Int32Array(this.pNumCells + 1);
    this.cellParticleIds = new Int32Array(maxParticles);
    this.numParticles = 0;
  }

  integrateParticles(dt, gravity) {
    for (let i = 0; i < this.numParticles; i++) {
      this.particleVel[2 * i + 1] += dt * gravity;
      this.particlePos[2 * i] += this.particleVel[2 * i] * dt;
      this.particlePos[2 * i + 1] += this.particleVel[2 * i + 1] * dt;
    }
  }

  pushParticlesApart(numIters) {
    const colorDiffusionCoeff = 0.001;
    this.numCellParticles.fill(0);
    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i],
          y = this.particlePos[2 * i + 1];
      let xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      let yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      let cellNr = xi * this.pNumY + yi;
      this.numCellParticles[cellNr]++;
    }
    let first = 0;
    for (let i = 0; i < this.pNumCells; i++) {
      first += this.numCellParticles[i];
      this.firstCellParticle[i] = first;
    }
    this.firstCellParticle[this.pNumCells] = first;
    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i],
          y = this.particlePos[2 * i + 1];
      let xi = clamp(Math.floor(x * this.pInvSpacing), 0, this.pNumX - 1);
      let yi = clamp(Math.floor(y * this.pInvSpacing), 0, this.pNumY - 1);
      let cellNr = xi * this.pNumY + yi;
      this.firstCellParticle[cellNr]--;
      this.cellParticleIds[this.firstCellParticle[cellNr]] = i;
    }
    const minDist = 2.0 * this.particleRadius;
    const minDist2 = minDist * minDist;
    for (let iter = 0; iter < numIters; iter++) {
      for (let i = 0; i < this.numParticles; i++) {
        let px = this.particlePos[2 * i],
            py = this.particlePos[2 * i + 1];
        let pxi = Math.floor(px * this.pInvSpacing);
        let pyi = Math.floor(py * this.pInvSpacing);
        let x0 = Math.max(pxi - 1, 0);
        let y0 = Math.max(pyi - 1, 0);
        let x1 = Math.min(pxi + 1, this.pNumX - 1);
        let y1 = Math.min(pyi + 1, this.pNumY - 1);
        for (let xi = x0; xi <= x1; xi++) {
          for (let yi = y0; yi <= y1; yi++) {
            let cellNr = xi * this.pNumY + yi;
            let first = this.firstCellParticle[cellNr],
                last = this.firstCellParticle[cellNr + 1];
            for (let j = first; j < last; j++) {
              let id = this.cellParticleIds[j];
              if (id === i) continue;
              let qx = this.particlePos[2 * id],
                  qy = this.particlePos[2 * id + 1];
              let dx = qx - px,
                  dy = qy - py;
              let d2 = dx * dx + dy * dy;
              if (d2 > minDist2 || d2 === 0.0) continue;
              let d = Math.sqrt(d2),
                  s = 0.5 * (minDist - d) / d;
              dx *= s; dy *= s;
              this.particlePos[2 * i] -= dx;
              this.particlePos[2 * i + 1] -= dy;
              this.particlePos[2 * id] += dx;
              this.particlePos[2 * id + 1] += dy;
              for (let k = 0; k < 3; k++) {
                let color0 = this.particleColor[3 * i + k],
                    color1 = this.particleColor[3 * id + k],
                    color = (color0 + color1) * 0.5;
                this.particleColor[3 * i + k] += (color - color0) * colorDiffusionCoeff;
                this.particleColor[3 * id + k] += (color - color1) * colorDiffusionCoeff;
              }
            }
          }
        }
      }
    }
  }

  handleParticleCollisions(obstacleX, obstacleY, obstacleRadius, obstacleVelX, obstacleVelY) {
    let h = 1.0 / this.fInvSpacing,
        r = this.particleRadius,
        minDist = obstacleRadius + r,
        minDist2 = minDist * minDist,
        minX = h + r, maxX = (this.fNumX - 1) * h - r,
        minY = h + r, maxY = (this.fNumY - 1) * h - r;
    for (let i = 0; i < this.numParticles; i++) {
      let x = this.particlePos[2 * i],
          y = this.particlePos[2 * i + 1];
      let dx = x - obstacleX, dy = y - obstacleY, d2 = dx * dx + dy * dy;
      if (d2 < minDist2) {
        this.particleVel[2 * i] = obstacleVelX;
        this.particleVel[2 * i + 1] = obstacleVelY;
      }
      if (x < minX) { x = minX; this.particleVel[2 * i] = 0.0; }
      if (x > maxX) { x = maxX; this.particleVel[2 * i] = 0.0; }
      if (y < minY) { y = minY; this.particleVel[2 * i + 1] = 0.0; }
      if (y > maxY) { y = maxY; this.particleVel[2 * i + 1] = 0.0; }
      this.particlePos[2 * i] = x;
      this.particlePos[2 * i + 1] = y;
    }
  }

  updateParticleDensity() {
    const n = this.fNumY, h = this.h, h1 = this.fInvSpacing, h2 = 0.5 * h;
    let d = this.particleDensity;
    d.fill(0.0);
    for (let i = 0; i < this.numParticles; i++) {
      let x = clamp(this.particlePos[2 * i], h, (this.fNumX - 1) * h);
      let y = clamp(this.particlePos[2 * i + 1], h, (this.fNumY - 1) * h);
      let x0 = Math.floor((x - h2) * h1), tx = ((x - h2) - x0 * h) * h1;
      let x1 = Math.min(x0 + 1, this.fNumX - 2);
      let y0 = Math.floor((y - h2) * h1), ty = ((y - h2) - y0 * h) * h1;
      let y1 = Math.min(y0 + 1, this.fNumY - 2);
      let sx = 1.0 - tx, sy = 1.0 - ty;
      d[x0 * n + y0] += sx * sy;
      d[x1 * n + y0] += tx * sy;
      d[x1 * n + y1] += tx * ty;
      d[x0 * n + y1] += sx * ty;
    }
    if (this.particleRestDensity === 0.0) {
      let sum = 0.0, numFluidCells = 0;
      for (let i = 0; i < this.fNumCells; i++) {
        if (this.cellType[i] === FLUID_CELL) {
          sum += d[i];
          numFluidCells++;
        }
      }
      if (numFluidCells > 0)
        this.particleRestDensity = sum / numFluidCells;
    }
  }

  transferVelocities(toGrid, flipRatio) {
    const n = this.fNumY, h = this.h, h1 = this.fInvSpacing, h2 = 0.5 * h;
    if (toGrid) {
      this.prevU.set(this.u);
      this.prevV.set(this.v);
      this.du.fill(0.0);
      this.dv.fill(0.0);
      this.u.fill(0.0);
      this.v.fill(0.0);
      for (let i = 0; i < this.fNumCells; i++)
        this.cellType[i] = this.s[i] === 0.0 ? SOLID_CELL : AIR_CELL;
      for (let i = 0; i < this.numParticles; i++) {
        let x = this.particlePos[2 * i],
            y = this.particlePos[2 * i + 1];
        let xi = clamp(Math.floor(x * h1), 0, this.fNumX - 1),
            yi = clamp(Math.floor(y * h1), 0, this.fNumY - 1);
        let cellNr = xi * n + yi;
        if (this.cellType[cellNr] === AIR_CELL)
          this.cellType[cellNr] = FLUID_CELL;
      }
    }
    for (let component = 0; component < 2; component++) {
      let dx = component === 0 ? 0.0 : h2;
      let dy = component === 0 ? h2 : 0.0;
      let f = component === 0 ? this.u : this.v;
      let prevF = component === 0 ? this.prevU : this.prevV;
      let dArr = component === 0 ? this.du : this.dv;
      for (let i = 0; i < this.numParticles; i++) {
        let x = clamp(this.particlePos[2 * i], h, (this.fNumX - 1) * h);
        let y = clamp(this.particlePos[2 * i + 1], h, (this.fNumY - 1) * h);
        let x0 = Math.min(Math.floor((x - dx) * h1), this.fNumX - 2);
        let tx = ((x - dx) - x0 * h) * h1;
        let x1 = Math.min(x0 + 1, this.fNumX - 2);
        let y0 = Math.min(Math.floor((y - dy) * h1), this.fNumY - 2);
        let ty = ((y - dy) - y0 * h) * h1;
        let y1 = Math.min(y0 + 1, this.fNumY - 2);
        let sx = 1.0 - tx, sy = 1.0 - ty;
        let d0 = sx * sy, d1 = tx * sy, d2 = tx * ty, d3 = sx * ty;
        let nr0 = x0 * n + y0, nr1 = x1 * n + y0, nr2 = x1 * n + y1, nr3 = x0 * n + y1;
        if (toGrid) {
          let pv = this.particleVel[2 * i + component];
          f[nr0] += pv * d0;  dArr[nr0] += d0;
          f[nr1] += pv * d1;  dArr[nr1] += d1;
          f[nr2] += pv * d2;  dArr[nr2] += d2;
          f[nr3] += pv * d3;  dArr[nr3] += d3;
        } else {
          let offset = component === 0 ? n : 1;
          let valid0 = (this.cellType[nr0] !== AIR_CELL || this.cellType[nr0 - offset] !== AIR_CELL) ? 1.0 : 0.0;
          let valid1 = (this.cellType[nr1] !== AIR_CELL || this.cellType[nr1 - offset] !== AIR_CELL) ? 1.0 : 0.0;
          let valid2 = (this.cellType[nr2] !== AIR_CELL || this.cellType[nr2 - offset] !== AIR_CELL) ? 1.0 : 0.0;
          let valid3 = (this.cellType[nr3] !== AIR_CELL || this.cellType[nr3 - offset] !== AIR_CELL) ? 1.0 : 0.0;
          let v = this.particleVel[2 * i + component];
          let d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;
          if (d > 0.0) {
            let picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] +
                        valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
            let corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) +
                        valid1 * d1 * (f[nr1] - prevF[nr1]) +
                        valid2 * d2 * (f[nr2] - prevF[nr2]) +
                        valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
            let flipV = v + corr;
            // Use 20% PIC and 80% FLIP:
            this.particleVel[2 * i + component] = 0.2 * picV + 0.8 * flipV;
          }
        }
      }
      if (toGrid) {
        for (let i = 0; i < f.length; i++) {
          if (dArr[i] > 0.0) f[i] /= dArr[i];
        }
        for (let i = 0; i < this.fNumX; i++) {
          for (let j = 0; j < this.fNumY; j++) {
            let idx = i * n + j;
            if (this.cellType[idx] === SOLID_CELL || (i > 0 && this.cellType[(i - 1) * n + j] === SOLID_CELL))
              this.u[idx] = this.prevU[idx];
            if (this.cellType[idx] === SOLID_CELL || (j > 0 && this.cellType[i * n + j - 1] === SOLID_CELL))
              this.v[idx] = this.prevV[idx];
          }
        }
      }
    }
  }

  // Fluid cells rendered as grayscale (higher density = darker).
  setSciColor(cellNr, val, minVal, maxVal) {
    let t = (val - minVal) / (maxVal - minVal);
    t = clamp(t, 0, 1);
    let gray = 1.0 - t;
    this.cellColor[3 * cellNr] = gray;
    this.cellColor[3 * cellNr + 1] = gray;
    this.cellColor[3 * cellNr + 2] = gray;
  }

  updateCellColors() {
    this.cellColor.fill(0.0);
    for (let i = 0; i < this.fNumCells; i++) {
      if (this.cellType[i] === SOLID_CELL) {
        // Draw obstacle (solid) in black.
        this.cellColor[3 * i] = 0.0;
        this.cellColor[3 * i + 1] = 0.0;
        this.cellColor[3 * i + 2] = 0.0;
      } else if (this.cellType[i] === FLUID_CELL) {
        let d = this.particleDensity[i];
        if (this.particleRestDensity > 0.0) d /= this.particleRestDensity;
        this.setSciColor(i, d, 1.5, 0.0);
      }
    }
  }

  solveIncompressibility(numIters, dt, overRelaxation, compensateDrift = true) {
    this.p.fill(0.0);
    this.prevU.set(this.u);
    this.prevV.set(this.v);
    const n = this.fNumY, cp = this.density * this.h / dt;
    for (let iter = 0; iter < numIters; iter++) {
      for (let i = 1; i < this.fNumX - 1; i++) {
        for (let j = 1; j < this.fNumY - 1; j++) {
          if (this.cellType[i * n + j] !== FLUID_CELL) continue;
          let center = i * n + j,
              left = (i - 1) * n + j,
              right = (i + 1) * n + j,
              bottom = i * n + j - 1,
              top = i * n + j + 1;
          let s = this.s[left] + this.s[right] + this.s[bottom] + this.s[top];
          if (s === 0.0) continue;
          let div = this.u[right] - this.u[center] + this.v[top] - this.v[center];
          if (this.particleRestDensity > 0.0 && compensateDrift) {
            let k = 1.0;
            let compression = this.particleDensity[i * n + j] - this.particleRestDensity;
            if (compression > 0.0) div -= k * compression;
          }
          let pVal = -div / s * overRelaxation;
          this.p[center] += cp * pVal;
          this.u[center] -= this.s[left] * pVal;
          this.u[right] += this.s[right] * pVal;
          this.v[center] -= this.s[bottom] * pVal;
          this.v[top] += this.s[top] * pVal;
        }
      }
    }
  }

  simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, overRelaxation, compensateDrift, separateParticles, obstacleX, obstacleY, obstacleRadius, obstacleVelX, obstacleVelY) {
    const numSubSteps = 1, sdt = dt / numSubSteps;
    for (let step = 0; step < numSubSteps; step++) {
      this.integrateParticles(sdt, gravity);
      if (separateParticles) this.pushParticlesApart(numParticleIters);
      this.handleParticleCollisions(obstacleX, obstacleY, obstacleRadius, obstacleVelX, obstacleVelY);
      this.transferVelocities(true);
      this.updateParticleDensity();
      this.solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
      this.transferVelocities(false, flipRatio);
    }
    this.updateCellColors();
  }
}