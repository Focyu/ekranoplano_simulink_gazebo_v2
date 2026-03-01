# Ekranoplano Simulink Gazebo V2 🚁🌊

Este repositorio contiene el entorno de simulación tridimensional y las arquitecturas de control PID desarrolladas para un Ekranoplano. El sistema integra el modelado dinámico en **MATLAB/Simulink**, la generación automática de código C++ mediante **Simulink Coder**, y la simulación física en **Gazebo** bajo el framework de **ROS 2**.

## 📁 Estructura del Repositorio

```text
ekranoplano_simulink_gazebo_v2/
├── MODELO GAZEBO/         # Modelos 3D (.sdf, .urdf) y mundos de simulación para Gazebo
├── ekranoplano_sim/       # Paquete de ROS (nodos, launch files, topics)
├── matlab_scripts/        # Scripts de inicialización (.m) y cálculos matemáticos (ej. rot_body_to_ned.m)
└── simulink_coder/        # Código C++ autogenerado por Simulink para los controladores (PID, etc.)

