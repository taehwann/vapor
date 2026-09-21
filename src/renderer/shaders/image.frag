#version 330 core
in vec2 uv;
out vec4 color;
uniform sampler2D u_Density;
uniform float u_Exposure;
uniform bool u_HistoricalPalette;
void main() {
    float density = max(texture(u_Density, uv).r, 0.0);
    if (u_HistoricalPalette) {
        color = vec4(mix(vec3(.03,.04,.07), vec3(.78,.80,.87), clamp(density*2.5,0.,1.)),1.);
        return;
    }
    float brightness = 1.0 - exp(-u_Exposure * density);
    color = vec4(brightness * vec3(0.8, 0.9, 1.0), 1.0);
}
