#' ODE Solution
#'
#' This function solves an ODE system
#' @param user created parameters
#' @return solution of an ODE function system
#' @export



# 定义模型
model <- function(t, y, param) {
  # 从参数列表中获得5种状态 SEIR
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  R <- y[4]
  I2 <- y[5] # 添加二次感染状态变量

  N <- param["N"]
  # 参数列表中提取模型参数
  beta <- param['beta']
  mu <- param['mu']
  gamma <- param['gamma']
  lamda <- param['lamda']
  delta <- param['delta'] # 新参数，表示二次感染率

  # 基于模型参数和状态信息创建传染并传播模型
  dSt <- mu * (N - S) - beta * S * (I1 + I2) / N # 易感者
  dEt <- beta * S * (I1 + I2) / N - mu * E - lamda * E # 暴露者（潜伏者）
  dI1t <- -(mu + delta ) * I1 + lamda * E  # 感染者
  dI2t <- -(mu + gamma) * I2 + delta * I1  # 二次感染者
  dRt <- - mu * R + gamma * I2 # 康复者

  # 汇总模型结果
  outcome <- c(dSt, dEt, dI1t, dRt, dI2t)
  list(outcome)
}

# 设置仿真参数
times <- seq(0, 200, by = 1 / 7)
param <- c(mu = 0.000, lamda = 0.03, beta = 4, gamma = 0.1, delta = 0.15, N = 1) # 添加 delta 参数
init <- c(S = 0.9999, E = 0.00000, I1 = 0.00002, R = 0, I2 = 0) # 初始值中加入 I2

# 微分方程求解函数
result <- deSolve::ode(y = init, times = times, func = model, parms = param)
result <- as.data.frame(result)

# 尾部数据变化情况
tail(round(result, 3.6), 10)

# 绘制图形
#' @export
seirplot <- ggplot2::ggplot(data = result) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = S, col = 'S'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I1, col = 'I1'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = R, col = 'R'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = E, col = 'E'), lwd = 2) +
  ggplot2::geom_line(ggplot2::aes(x = time, y = I2, col = 'I2'), lwd = 2) + # 添加二次感染的曲线
  ggplot2::labs(x = "时间", y = "比率") +
  ggplot2::scale_color_manual(name = "传染病模型仿真", values = c("S" = "orange", "E" = "purple", "I1" = "pink", "R" = "skyblue", "I2" = "green")) # 修改颜色映射

# 绘制图形对象
seirplot


# 保存为矢量图
ggplot2::ggsave(seirplot, file = "file3.pdf", width=7, height=6)
ggplot2::ggsave(seirplot, file = "file3.svg",  width=7, height=6)
