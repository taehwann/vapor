#version 330 core
out vec2 uv;
void main() {
    // Oversized triangle covers the viewport without a diagonal seam.
    vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    uv = p;
    gl_Position = vec4(2.0 * p - 1.0, 0.0, 1.0);
}
