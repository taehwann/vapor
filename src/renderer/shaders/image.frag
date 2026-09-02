#version 330 core
in vec2 uv;
out vec4 color;
uniform sampler2D u_Density;
uniform float u_Exposure;
void main() {
    float density = max(texture(u_Density, uv).r, 0.0);
    float brightness = 1.0 - exp(-u_Exposure * density);
    color = vec4(brightness * vec3(0.8, 0.9, 1.0), 1.0);
}
