# =============================================================================
# 对比 Theis 解析解 与 2D 瞬态有限差分数值解
#
# 目的：验证 Ftransient2dsim 在均质承压含水层中是否能正确逼近 Theis 解析解。
#       如果二者偏差较大，说明 2D 瞬态流代码可能存在 bug。
#
# 设计思路：
#   1. 用较大的区域 (domain) + 较细的网格，使数值解逼近"无限延伸"含水层
#   2. 观测孔设置在离抽水井较近的位置（避免边界效应影响）
#   3. 时间跨度覆盖从早期到接近稳态
#   4. 对比降深曲线 (semilog 和 log-log)
# =============================================================================

library(HydroTomo)
library(ggplot2)

# ---------------------------------------------------------------------------
# 含水层参数 (与 Theis 解析解一致)
# ---------------------------------------------------------------------------
Q  <- 1.3e-3       # 抽水流量 [m³/s]
TT <- 1.5e-3       # 导水系数 [m²/s]
SS <- 2.0e-5       # 储水系数 [-]
rw <- 200          # 观测孔距离 [m]

# ---------------------------------------------------------------------------
# 1. Theis 解析解 (用于绘制光滑参考曲线)
# ---------------------------------------------------------------------------
t_theis <- 10^seq(1, 5, length.out = 200)
s_theis <- theis_drawdown(Q = Q, r = rw, T = TT, S = SS, t = t_theis)
s_cj    <- cooper_jacob_drawdown(Q, rw, TT, SS, t_theis)

# ---------------------------------------------------------------------------
# 2. 2D 瞬流数值解
# ---------------------------------------------------------------------------
# 区域要足够大，使得在观测时间范围内边界效应不明显
# 用 L²/(4*T*t_max) 来估计影响半径 ~ sqrt(4*T*t_max / S)
# 取 t_max = 10^5 s, Ri ~ sqrt(4*1.5e-3*1e5/2e-5) ≈ sqrt(3e4) ≈ 173 m
# 所以区域至少 2*Ri ≈ 350 m 宽才能覆盖影响范围
# 我们用 800×800 m，抽水井放中心，观测井距中心 200 m
nx <- 80
ny <- 80
domain <- c(nx, ny, 0, 800, 0, 800)
grid <- GenGrid(domain)

# 抽水井在区域中心
xp <- 400; yp <- 400
Qinf <- data.frame(Qp = Q, x = xp, y = yp)

# 时间向量：seq(0, 1000, by = 1)，密集的时间步长减少插值误差，对数坐标下更光滑
times <- seq(0, 1000, by = 1)

# 观测时刻：从数值解中提取
obs_time <- c(1, 3, 5, 10, 30, 50, 100, 300, 500,
              1000, 3000, 5000, 10000, 30000, 50000, 100000)
obs_time <- obs_time[obs_time <= max(times)]  # 只保留在模拟时间范围内的

cat(sprintf("网格: %d x %d, 区域: 800 x 800 m², dx = dy = %.1f m\n", nx, ny, grid$dx))
cat(sprintf("抽水井: (%.0f, %.0f), 观测孔: (%.0f, %.0f) r = %.0f m\n",
            xp, yp, xp + rw, yp, rw))
cat(sprintf("时间: %d 步 (0 到 %.0f s, dt = 1 s)\n", length(times), max(times)))
cat("正在运行 2D 瞬态模拟...\n")

res <- Ftransient2dsim(domain = domain, grid = grid,
                        TT = TT, SS = SS, Qinf = Qinf, times = times,
                        lrw = 2000000)

cat("2D 模拟完成。\n")

# 提取观测孔处的降深
Oinf <- data.frame(data = NA, x = xp + rw, y = yp, time = obs_time)
Oinf <- samDataTr(Oinf = Oinf, grid = res$grid, result_tr = res)
s_2d <- -Oinf$data  # 转为降深（正值）

# 对应时刻的 Theis 解
s_theis_obs <- theis_drawdown(Q = Q, r = rw, T = TT, S = SS, t = obs_time)

# ---------------------------------------------------------------------------
# 3. 构建 ggplot2 数据
# ---------------------------------------------------------------------------
df_theis <- data.frame(t = t_theis, s = s_theis)
df_cj    <- data.frame(t = t_theis, s = s_cj)
df_2d    <- data.frame(t = obs_time, s = s_2d)
df_err   <- data.frame(
  t       = obs_time,
  rel_err = (s_2d - s_theis_obs) / s_theis_obs * 100
)

# 降深场数据 (t = max time)
h_last <- matrix(res$out[nrow(res$out), -1], nx, ny)
df_field <- expand.grid(x = grid$xmid, y = grid$ymid)
df_field$s <- -as.vector(h_last)

# 井位标注
df_wells <- data.frame(
  x = c(xp, xp + rw),
  y = c(yp, yp),
  label = c("Pumping well", "Observation well")
)

# ---------------------------------------------------------------------------
# 4. 绘图
# ---------------------------------------------------------------------------

# (a) Semilog 降深曲线
p1 <- ggplot() +
  geom_line(data = df_theis, aes(t, s, color = "Theis analytical"), linewidth = 1) +
  geom_line(data = df_cj, aes(t, s, color = "Cooper-Jacob approx"),
            linetype = "dashed", linewidth = 0.6) +
  geom_point(data = df_2d, aes(t, s, color = "2D numerical"), size = 2.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c("Theis analytical" = "steelblue",
               "Cooper-Jacob approx" = "gray50",
               "2D numerical" = "red")
  ) +
  labs(x = "Time t (s)", y = "Drawdown s (m)",
       title = "(a) Semilog drawdown curve", color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "gray80"))

# (b) Log-log 降深曲线
p2 <- ggplot() +
  geom_line(data = df_theis, aes(t, s, color = "Theis analytical"), linewidth = 1) +
  geom_point(data = df_2d, aes(t, s, color = "2D numerical"), size = 2.5) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(
    values = c("Theis analytical" = "steelblue", "2D numerical" = "red")
  ) +
  labs(x = "Time t (s)", y = "Drawdown s (m)",
       title = "(b) Log-log drawdown curve", color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "gray80"))

# (c) 相对偏差
p3 <- ggplot(df_err, aes(t, rel_err)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.6, color = "darkred") +
  geom_point(size = 2.5, color = "darkred") +
  scale_x_log10() +
  labs(x = "Time t (s)", y = "Relative error (%)",
       title = "(c) 2D numerical vs Theis relative error") +
  theme_bw(base_size = 13)

# (d) 降深场快照
p4 <- ggplot(df_field, aes(x, y, fill = s)) +
  geom_raster() +
  scale_fill_gradientn(colors = terrain.colors(64)) +
  geom_point(data = df_wells, aes(x, y, shape = label),
             size = 3, color = "black", fill = "white", inherit.aes = FALSE) +
  scale_shape_manual(values = c("Pumping well" = 8, "Observation well" = 1)) +
  coord_fixed(ratio = 1) +
  labs(x = "x (m)", y = "y (m)", fill = "s (m)",
       shape = NULL,
       title = sprintf("(d) 2D drawdown field  t = %.0f s", max(times))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

# 组合输出
ggsave("compare_theis_2d.pdf",
       gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2),
       width = 12, height = 10)

cat(sprintf("\n图片已保存: compare_theis_2d.pdf\n"))

# ---------------------------------------------------------------------------
# 统计输出
# ---------------------------------------------------------------------------
cat(sprintf("\n相对偏差统计:\n"))
cat(sprintf("  平均偏差: %.2f %%\n", mean(df_err$rel_err, na.rm = TRUE)))
cat(sprintf("  最大偏差: %.2f %%\n", max(abs(df_err$rel_err), na.rm = TRUE)))
cat(sprintf("  RMSE: %.4f m\n", sqrt(mean((s_2d - s_theis_obs)^2, na.rm = TRUE))))

# ---------------------------------------------------------------------------
# 5. 检查 u < 0.03 范围内的拟合
# ---------------------------------------------------------------------------
u_obs <- (SS * rw^2) / (4 * TT * obs_time)
idx_cj <- u_obs < 0.03
if (any(idx_cj)) {
  s_cj_obs <- cooper_jacob_drawdown(Q, rw, TT, SS, obs_time[idx_cj])
  cat(sprintf("\n在 u < 0.03 范围内 (%d 个点):\n", sum(idx_cj)))
  cat(sprintf("  2D vs Theis  RMSE: %.4f m\n",
      sqrt(mean((s_2d[idx_cj] - s_theis_obs[idx_cj])^2))))
  cat(sprintf("  2D vs CJ     RMSE: %.4f m\n",
      sqrt(mean((s_2d[idx_cj] - s_cj_obs)^2))))
}

# ---------------------------------------------------------------------------
# 6. 数值解收敛性检查 —— 不同网格精度对比
# ---------------------------------------------------------------------------
cat("\n--- 网格收敛性检查 ---\n")
grid_sizes <- c(40, 80, 120)
results <- list()
for (nx_test in grid_sizes) {
  ny_test <- nx_test
  dom_test <- c(nx_test, ny_test, 0, 800, 0, 800)
  grid_test <- GenGrid(dom_test)
  times_test <- seq(0, 1000, by = 1)
  res_test <- Ftransient2dsim(domain = dom_test, grid = grid_test,
                               TT = TT, SS = SS, Qinf = Qinf,
                               times = times_test, lrw = 2000000)
  obs_test <- c(10, 100, 1000)
  Oinf_test <- data.frame(data = NA, x = xp + rw, y = yp, time = obs_test)
  Oinf_test <- samDataTr(Oinf = Oinf_test, grid = res_test$grid,
                          result_tr = res_test)
  s_2d_test <- -Oinf_test$data
  s_ref <- theis_drawdown(Q, rw, TT, SS, obs_test)
  rmse_test <- sqrt(mean((s_2d_test - s_ref)^2))
  results[[as.character(nx_test)]] <- rmse_test
  cat(sprintf("  网格 %d×%d: RMSE = %.4f m\n", nx_test, ny_test, rmse_test))
}
