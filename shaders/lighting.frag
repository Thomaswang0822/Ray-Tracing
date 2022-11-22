#version 330 core

in vec4 position; // raw position in the model coord
in vec3 normal;   // raw normal in the model coord

uniform mat4 modelview; // from model coord to eye coord
uniform mat4 view;      // from world coord to eye coord

// Material parameters
uniform vec4 ambient;
uniform vec4 diffuse;
uniform vec4 specular;
uniform vec4 emision;
uniform float shininess;

// Light source parameters
const int maximal_allowed_lights = 10;
uniform bool enablelighting;
uniform int nlights;
uniform vec4 lightpositions[ maximal_allowed_lights ]; // in world coord
uniform vec4 lightcolors[ maximal_allowed_lights ];

// Output the frag color
out vec4 fragColor;


void main (void){
    if (!enablelighting){
        // Default normal coloring (you don't need to modify anything here)
        vec3 N = normalize(normal);
        fragColor = vec4(0.5f*N + 0.5f , 1.0f);
    } else {
        
        // HW3: You will compute the lighting here.
        fragColor = emision;
        
        // Find v
        // We don't know eye position in model coor
        // So we turn object point position to camera matrix coor
        vec3 eye = vec3(0.0f, 0.0f, 0.0f);
        vec4 p_cam = modelview * position;
        vec3 v = normalize(eye - vec3(p_cam)/p_cam[3]);

        // Find n in model coor to camera coor by 3x3 VM ^ (-T)
        mat3 vm3 = mat3(modelview); 
        vm3 = transpose(inverse(vm3));
        vec3 n = normalize(vm3 * normal);
        for (int j=0; j<nlights; j++) {
            vec4 inside = ambient;
            // Find lj, unit vector to the light
            vec4 light_cam = view * lightpositions[j];
            vec3 lj = normalize(
                vec3(light_cam) * p_cam[3] - 
                vec3(p_cam) * light_cam[3]
            );
            if (dot(n, lj) > 0) {
                inside = inside + diffuse * dot(n, lj);
            }
            // find hj
            vec3 hj = normalize(v + lj);
            if (dot(n, hj) > 0) {
                inside = inside + specular * pow( dot(n, hj), shininess); 
            }
            inside = inside * lightcolors[j];

            fragColor = fragColor + inside;
        }
        
    }
}
