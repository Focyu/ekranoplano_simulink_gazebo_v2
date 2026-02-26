#!/bin/bash

# Imprimir un mensaje de inicio
echo "🚀 Iniciando simulación del Ekranoplano en Gazebo (ROS 2 Humble)..."

# 1. Cargar el entorno de ROS 2 (por si la terminal no lo hizo automáticamente)
source /opt/ros/humble/setup.bash

# 2. Entrar a la carpeta del proyecto
cd ~/Documentos/ekranoplano_sim

# 3. Decirle a Gazebo dónde está nuestro modelo 3D local
export GAZEBO_MODEL_PATH=${PWD}/models:$GAZEBO_MODEL_PATH

# 4. Lanzar Gazebo con el puente de ROS y nuestro mundo
ros2 launch gazebo_ros gazebo.launch.py world:=vuelo.world
