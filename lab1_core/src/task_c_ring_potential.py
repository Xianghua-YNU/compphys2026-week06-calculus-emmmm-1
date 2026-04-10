import numpy as np

def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    """
    用离散积分计算单点电势
    """
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    d_phi = 2 * np.pi / n_phi

    r2 = (x - a * np.cos(phi)) ** 2 + (y - a * np.sin(phi)) ** 2 + z ** 2
    integrand = 1.0 / np.sqrt(r2 + 1e-12)

    # 复合梯形积分
    integral = d_phi * (np.sum(integrand) - 0.5 * integrand[0] - 0.5 * integrand[-1])
    return q / (2 * np.pi) * integral

def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    """
    ✅ 修复索引错误：兼容一维数组 + 二维网格
    自动适配测试用的一维输入 和 绘图用的二维网格
    """
    # 保存原始形状，展平为一维计算
    shape = y_grid.shape
    y_flat = y_grid.ravel()
    z_flat = z_grid.ravel()
    
    V_flat = np.zeros_like(y_flat, dtype=float)
    for i in range(len(y_flat)):
        V_flat[i] = ring_potential_point(x0, y_flat[i], z_flat[i], a, q, n_phi)
    
    # 恢复原始形状
    return V_flat.reshape(shape)

def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    return q / np.sqrt(a * a + z * z)
