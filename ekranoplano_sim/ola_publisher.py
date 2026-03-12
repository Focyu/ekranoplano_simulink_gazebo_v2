#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64
import math
import time

class OlaPublisher(Node):
    def __init__(self):
        super().__init__('ola_publisher')
        self.pub = self.create_publisher(Float64, '/ocean/wave_height', 10)
        self.timer = self.create_timer(0.05, self.tick)
        self.A     = 0.20
        self.T     = 6.5
        self.omega = 2 * math.pi / self.T
        self.t0    = time.time()
        self.get_logger().info('🌊 Ola publisher iniciado — Sea State 2 (A=0.20m, T=6.5s)')

    def tick(self):
        t   = time.time() - self.t0
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

