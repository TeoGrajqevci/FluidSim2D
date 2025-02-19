// common.js
export const canvas = document.getElementById("myCanvas");
export const gl = canvas.getContext("webgl");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
canvas.focus();

export const simHeight = 3.0;
export const cScale = canvas.height / simHeight;
export const simWidth = canvas.width / cScale;

export const FLUID_CELL = 0, AIR_CELL = 1, SOLID_CELL = 2;

export function clamp(x, min, max) {
  return Math.max(min, Math.min(max, x));
}
