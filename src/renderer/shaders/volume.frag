#version 330 core
            in vec3 vWorldPos;
            out vec4 FragColor;
            uniform sampler3D u_Volume;
            uniform float u_StepScale;
            uniform float u_BoxSize;
            uniform float u_GridRes;
            uniform float u_AlphaMul;
            uniform vec3 u_CamPos;
            uniform vec3 u_LightDir;
            uniform float u_ShadowStr;
            uniform float u_ShadowStep;
            void main() {
                vec3 ro = u_CamPos;
                vec3 rd = normalize(vWorldPos - ro);
                vec3 t0 = (vec3(0.0) - ro) / rd;
                vec3 t1 = (vec3(u_BoxSize) - ro) / rd;
                vec3 tn = min(t0, t1), tf = max(t0, t1);
                float tnear = max(max(tn.x, tn.y), tn.z);
                float tfar  = min(min(tf.x, tf.y), tf.z);
                if (tnear < 0.0) tnear = 0.0;
                if (tnear >= tfar) discard;
                float step = u_StepScale / u_GridRes;
                vec3 bg = vec3(0.03, 0.04, 0.07);
                vec3 smokeCol = vec3(0.9, 0.85, 0.75);
                vec4 col = vec4(0.0);
                bool hitSmoke = false;
                for (float t = tnear; t < tfar; t += step) {
                    float d = texture(u_Volume, (ro + rd * t) / u_BoxSize).r;
                    if (d > 0.001) {
                        hitSmoke = true;
                        float shadow = 1.0;
                        if (u_ShadowStr > 0.0) {
                            vec3 spos = ro + rd * t;
                            float sa = 0.0;
                            for (float st = u_ShadowStep; st < u_BoxSize * 2.0; st += u_ShadowStep) {
                                vec3 sp = spos + u_LightDir * st;
                                if (sp.x < 0.0 || sp.x > u_BoxSize || sp.y < 0.0 || sp.y > u_BoxSize || sp.z < 0.0 || sp.z > u_BoxSize) break;
                                sa += texture(u_Volume, sp / u_BoxSize).r * u_ShadowStep * u_AlphaMul * u_ShadowStr;
                                if (sa > 3.0) break;
                            }
                            shadow = exp(-sa);
                        }
                        float alpha = clamp(d * step * u_AlphaMul, 0.0, 1.0);
                        col.rgb += (1.0 - col.a) * smokeCol * alpha * shadow;
                        col.a += (1.0 - col.a) * alpha;
                    }
                    if (col.a > 0.99) break;
                }
                if (hitSmoke) {
                    FragColor = vec4(bg * (1.0 - col.a) + col.rgb, 1.0);
                } else {
                    discard;
                }
            }
