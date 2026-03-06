# Ekranoplano Simulink Gazebo V2 🚁🌊

Este repositorio contiene el entorno de simulación tridimensional y las arquitecturas de control PID desarrolladas para un Ekranoplano. El sistema integra el modelado dinámico en **MATLAB/Simulink**, la generación automática de código C++ mediante **Simulink Coder**, y la simulación física en **Gazebo** bajo el framework de **ROS 2**.

## 📁 Estructura del Repositorio

```text
ekranoplano_simulink_gazebo_v2/
├── MODELO GAZEBO/               # Modelos 3D (.sdf, .urdf) y mundos para Gazebo
├── ekranoplano_sim/             # Paquete principal de ROS 2 (launch, mundos, config)
├── matlab_scripts/              # Scripts matemáticos e inicialización de ganancias
└── pid_controllers_simulink/    # Controladores C++ autogenerados por Simulink
    └── pid_control_v1/          # Paquete ROS 2 del nodo controlador PID

