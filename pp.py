# 该代码实现了 L8-BO 与 Toluene 和 o-Xylene 的 UNIFAC 活度系数计算，并拟合 g12 曲线。
import numpy as np
from chemicals import unifac
from scipy.optimize import curve_fit

# 1. 定义摩尔体积 (单位: cm^3/mol, 来源参考前文查询结果)
V1_L8 = 1130.0  # L8-BO
V2_Tol = 106.3  # Toluene
V2_OXY = 120.6  # o-Xylene

# 2. 定义 UNIFAC 基团 (字典格式: {基团编号: 数量})
# 这里的基团编号遵循标准 UNIFAC (Hansen/Fredenslund)
groups_Tol = {9: 5, 10: 1}  # Toluene: 5*AC, 1*ACCH3
groups_OXY = {9: 4, 10: 2}  # o-Xylene: 4*AC, 2*ACCH3

# L8 的近似分解 (由于分子巨大且含复杂稠环，这里给出基于结构的估计)
# 实际计算中，建议根据分子式严格拆分，这里为演示流程
groups_L8 = {
    1: 6,   # CH3
    2: 40,  # CH2
    3: 4,   # CH
    9: 10,  # AC (芳香环)
    25: 4,  # C#N (氰基)
    18: 2,  # CH3CO (羰基简化)
    # 注意：标准UNIFAC可能缺乏某些稠环/硫杂环基团，通常用相似基团替代
}

def calculate_g12_profile(groups1, groups2, V1, V2, T=298.15):
    """
    计算 u2 从 0.01 到 0.99 对应的 g12 数组
    u2: 组分2(溶剂)的体积分数
    """
    u2_array = np.linspace(0.01, 0.99, 99)
    g12_array = []

    for u2 in u2_array:
        u1 = 1.0 - u2
        
        # --- 核心转换：体积分数(u)转为摩尔分数(x) ---
        # x1 = (u1/V1) / (u1/V1 + u2/V2)
        denom = (u1 / V1) + (u2 / V2)
        x1 = (u1 / V1) / denom
        x2 = (u2 / V2) / denom
        
        # --- 使用 UNIFAC 计算活度系数 gamma ---
        try:
            # unifac.gamma 返回每个组分的活度系数
            gammas = unifac.UF_gamma(T, [x1, x2], [groups1, groups2])
            gamma1, gamma2 = gammas[0], gammas[1]
            
            # --- 计算 g12 ---
            # 注意：u1, u2 是体积分数。根据用户公式：
            # g12 = (u1*ln(gamma1) + u2*ln(gamma2)) / (u1*u2)
            ln_gamma1 = np.log(gamma1)
            ln_gamma2 = np.log(gamma2)
            g12 = (u1 * ln_gamma1 + u2 * ln_gamma2) / (u1 * u2)
            g12_array.append(g12)
        except:
            # 异常处理：某些基团组合可能无参数
            g12_array.append(np.nan)
            
    return u2_array, np.array(g12_array)

def fit_g12(u2_array, g12_array):
    """
    四阶多项式拟合: p1*u2^4 + p2*u2^3 + p3*u2^2 + p4*u2 + p5
    """
    # 剔除无效值
    mask = ~np.isnan(g12_array)
    p = np.polyfit(u2_array[mask], g12_array[mask], 4)
    return p # 返回顺序为 [p1, p2, p3, p4, p5]

# --- 运行计算 ---
T_test = 298.15

# 情况 A: L8 (1) + Toluene (2)
u2_A, g12_A = calculate_g12_profile(groups_L8, groups_Tol, V1_L8, V2_Tol, T=T_test)
p_A = fit_g12(u2_A, g12_A)

# 情况 B: L8 (1) + o-Xylene (2)
u2_B, g12_B = calculate_g12_profile(groups_L8, groups_OXY, V1_L8, V2_OXY, T=T_test)
p_B = fit_g12(u2_B, g12_B)

# --- 最终输出 ---
print("L8 : Toluene 拟合系数 [p1, p2, p3, p4, p5]:")
print(p_A)
print("\nL8 : o-Xylene 拟合系数 [p1, p2, p3, p4, p5]:")
print(p_B)

# 验证公式输出示例
def g12_func(u2, p):
    return p[0]*u2**4 + p[1]*u2**3 + p[2]*u2**2 + p[3]*u2 + p[4]

print(f"\n验证: 当 u2=0.5 时, L8:Tol 的 g12 = {g12_func(0.5, p_A):.4f}")