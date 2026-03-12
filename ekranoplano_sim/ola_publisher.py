#!/usr/bin/env python3
"""
ola_publisher.py — Publica la altura de la ola del mar en tiempo real.
Sincroniza visualmente Gazebo con las olas calculadas en la MATLAB Function.
Mismos parámetros: A=0.20m, T=6.5s, lambda=50m (Sea State 2)
"""
import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64
import math
import time

class OlaPublisher(Node):
    def __init__(self):
        super().__init__('ola_publisher')

        self.pub = self.create_publisher(Float64, '/ocean/wave_height', 10)
        self.timer = self.create_timer(0.05, self.tick)  # 20 Hz

        # Parámetros — DEBEN coincidir exactamente con MATLAB Function
        self.A      = 0.20        # Amplitud [m]
        self.T      = 6.5         # Período [s]
        self.omega  = 2 * math.pi / self.T
        self.k      = 2 * math.pi / 50.0  # k = 2pi/lambda

        self.t0 = time.time()
        self.get_logger().info('🌊 Ola publisher iniciado — Sea State 2 (A=0.20m, T=6.5s)')

    def tick(self):
        t   = time.time() - self.t0
        # x10 ≈ posición NED del ekranoplano — simplificado en x=0 para la visual
        # La física exacta con x10 real la calcula Simulink internamente
        eta = self.A * math.cos(self.omega * t)

        msg = Float64()
        msg.data = float(eta)
        self.pub.publish(msg)

def main():
    rclpy.init()
    node = OlaPublisher()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()
