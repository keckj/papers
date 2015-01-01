
typedef struct Vec3 {
    float x;
    float y;
    float z;
} Vec3_s;

Vec3_s f(Vec3_s u) {
    Vec3_s v;
    v.x = 10.0f*(u.y - u.x);
    v.y = 28.0f*u.x - u.y - u.x*u.z;
    v.z = u.x*u.y - 8.0f/3.0f*u.z;
    return v;
}

Vec3_s RK4_step(Vec3_s u_k, float dt) {
    Vec3_s k1, k2, k3, k4;
    Vec3_s y2, y3, y4;
    Vec3_s v_k;

    k1 = f(u_k);
    y2.x = u_k.x + k1.x * dt/2.0f; 
    y2.y = u_k.y + k1.y * dt/2.0f; 
    y2.z = u_k.z + k1.z * dt/2.0f; 

    k2 = f(y2);
    y3.x = u_k.x + k2.x * dt/2.0f; 
    y3.y = u_k.y + k2.y * dt/2.0f; 
    y3.z = u_k.z + k2.z * dt/2.0f; 
    
    k3 = f(y3);
    y4.x = u_k.x + k3.x * dt; 
    y4.y = u_k.y + k3.y * dt; 
    y4.z = u_k.z + k3.z * dt; 

    k4 = f(y4);

    v_k.x = u_k.x + dt/6.0f*(k1.x + 2*k2.x + 2*k3.x + k4.x);
    v_k.y = u_k.y + dt/6.0f*(k1.y + 2*k2.y + 2*k3.y + k4.y);
    v_k.z = u_k.z + dt/6.0f*(k1.z + 2*k2.z + 2*k3.z + k4.z);
    
    return v_k;
}

