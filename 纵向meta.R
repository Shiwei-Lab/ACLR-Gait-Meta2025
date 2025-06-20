# 安装必要包（如果还没安装）
# install.packages(c("metafor", "clubSandwich", "splines"))
# install.packages("rms")      # 第一次才需要安装
library(metafor)
library(clubSandwich)
library(splines)
library(writexl)
library(rms) 
library(dplyr)
library(tidyverse)
# 1. 计算效应量及其方差（SMD = 标准化均值差）
dat <- read_xlsx("屈曲力矩/屈曲力矩数据集.xlsx")  # 替换为你自己的文件名和路径

effect_sizes <- escalc(
  measure = "SMD", 
  m1i = M_post_exp, sd1i = SD_post_exp, n1i = n_post_exp,   # 实验组
  m2i = M_post_ctrl, sd2i = SD_post_ctrl, n2i = n_post_ctrl, # 对照组
  data = dat
)
write_xlsx(effect_sizes, "屈曲力矩/屈曲力矩数据集.xlsx")

# 构建 CAR 结构的协方差矩阵，考虑不同时间点之间的自相关性
set.seed(123)
V <- impute_covariance_matrix(
  vi = dat$vi,
  cluster = dat$study_id,
  ti = dat$time_point,
  ar1 = 0.85,  # 假设时间点间自相关
  check_PD = TRUE,
  smooth_vi = TRUE,
  return_list = FALSE
)

# 空模型（仅随机效应，不含 moderators）
model_empty <- rma.mv(
  yi, V,
  random = ~ time_point | study_id,
  struct = "CAR",
  data = dat
)

# 3. 查看结果摘要
summary(model_empty)

#绘制森林图
# 添加研究标签（slab）
effect_sizes$slab <- paste0("", effect_sizes$study_id, " (", effect_sizes$time_point, " mo)")
forest(model_empty,
       slab = effect_sizes$slab,
       xlab = "Effect Size (SMD)",
       alim = c(-2.5, 2.5),
       cex = 0.9,         # 字体大小
       pch = 9,          # 方形点
       col = "pink",  # 线和点颜色
       lwd = 2.2,           # 置信区间线条粗细
       mlab = "RE Model", # 模型标签
       efac = 2      # 点的大小
)


# 自定义伪R²（pseudo-R²）用于衡量 moderator 模型比空模型解释的异质性比例）
r2_model <- function(model){
  null_model <- suppressWarnings(update(model, mods = NULL))
  pseudo_r2 <- 100 * (null_model$tau2 - model$tau2) / null_model$tau2
  return(pseudo_r2)
}
# 含不同 moderators 的 4 个模型
model_linear <- update(model_empty, mods = ~ time_point)
summary(model_linear)
model_log <- update(model_empty, mods = ~ log(time_point))
summary(model_log)
model_rcs3 <- update(model_empty, mods = ~ rcs(time_point, 3))
summary(model_rcs3)
model_rcs4 <- update(model_empty, mods = ~ rcs(time_point, 4)) 
summary(model_rcs4)

# 比较不同模型的 LogLik, AIC, BIC, tau², 伪R² 等指标
compare_models <- function(...){
  models <- list(...)
  names(models) <- c("Linear", "Log", "3 knot RCS", "4 knot RCS")

  comparison <- lapply(models, function(m){
    null <- suppressWarnings(update(m, mods = NULL))
    data.frame(
      LogLik = logLik(m),
      AIC = AIC(m),
      BIC = BIC(m),
      tau2 = m$tau2,
      pseudoR2 = 100 * (null$tau2 - m$tau2) / null$tau2
    )
  })

  do.call(rbind, comparison)
}
# 执行模型比较
compare_models(model_linear, model_log, model_rcs3, model_rcs4)

# 应用 clubSandwich 的 cluster-robust 修正（适用于小样本推断）
robust_rcs3 <- robust(model_rcs3, cluster = dat$study_id, clubSandwich = TRUE)
robust_linear <- robust(model_linear, cluster = dat$study_id, clubSandwich = TRUE)
robust_log <- robust(model_log, cluster = dat$study_id, clubSandwich = TRUE)
robust_rcs4 <- robust(model_rcs4, cluster = dat$study_id, clubSandwich = TRUE)
# 查看修正后的摘要
summary(robust_log)
summary(robust_linear)
summary(robust_rcs3)
summary(robust_rcs4)


#模型可视化
# 可视化 log 模型
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))
pred <- predict(model_log, newmods = log(new_data$time_point))

new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 4️⃣ 计算总样本数、研究数
total_k <- model_log$k
total_studies <- dat %>% filter(!is.na(vi)) %>% summarise(n = n_distinct(study_id)) %>% pull(n)
total_n <- dat %>%
  filter(!is.na(vi)) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = sum(n_total, na.rm = TRUE)) %>%
  pull(n)

# 5️⃣ 绘图
ggplot() +
  # 基准线
  geom_hline(yintercept = 0, colour = "dark grey") +
  
  # 每个 study_id 的线条
  geom_line(data = dat, aes(x = time_point, y = yi, group = study_id),
            alpha = 0.8, colour = "grey") +
  
  # 气泡散点
  geom_point(data = dat, aes(x = time_point, y = yi, size = n_post_exp),
             alpha = 0.3, colour = "black", fill = "grey") +
  
  # 外层更宽阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 0.35,
                                   ymax = pred + 0.35),
              fill = "black", alpha = 0.1) +
  
  # 内层窄阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 0.15,
                                   ymax = pred + 0.15),
              fill = "red", alpha = 0.2) +
  
  # 主曲线
  geom_line(data = new_data, aes(x = time_point, y = pred),
            colour = "red", linewidth = 1.3) +
  
  # x 轴（0–60）
  scale_x_continuous(
    breaks = seq(0,60, by = 6),
    labels = seq(0,60, by = 6),
    name = "Time since surgery (months)"
  ) +
  
  # y 轴（效应量）
  scale_y_continuous(name = "Effect size (yi)", limits = c(-3, 1)) +
  
  # 坐标范围
  coord_cartesian(xlim = c(0, 48)) +
  
  # 气泡大小
  scale_size(range = c(1, 6), name = "Participants (n)") +
  
  # 主题
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "plain"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_line(linewidth = rel(0.5), linetype = 2)
  ) +
  
  # 👇 注释
  annotate("text", x = 45, y = -2.6,
           label = paste0("K = ", total_k, " (", total_studies, " studies)"),
           hjust = 1, size = 3.3) +
  annotate("text", x = 40, y = -2.8,
           label = paste0("N = ", total_n),
           hjust = 1, size = 3.3)

# 可视化 linear 模型
# 新时间点数据
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 进行逐点预测
pred <- predict(model_linear, newmods = new_data$time_point)

# 将预测结果赋值到 new_data
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 4️⃣ 计算总样本数、研究数
total_k <- model_linear$k
total_studies <- dat %>% filter(!is.na(vi)) %>% summarise(n = n_distinct(study_id)) %>% pull(n)
total_n <- dat %>%
  filter(!is.na(vi)) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = sum(n_total, na.rm = TRUE)) %>%
  pull(n)

# 5️⃣ 绘图
ggplot() +
  # 基准线
  geom_hline(yintercept = 0, colour = "dark grey") +
  
  # 每个 study_id 的线条
  geom_line(data = dat, aes(x = time_point, y = yi, group = study_id),
            alpha = 0.8, colour = "grey") +
  
  # 气泡散点
  geom_point(data = dat, aes(x = time_point, y = yi, size = n_post_exp),
             alpha = 0.3, colour = "black", fill = "grey") +
  
  # 外层更宽阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 0.5,
                                   ymax = pred + 0.5),
              fill = "black", alpha = 0.1) +
  
  # 内层窄阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 0.15,
                                   ymax = pred + 0.15),
              fill = "red", alpha = 0.2) +
  
  # 主曲线
  geom_line(data = new_data, aes(x = time_point, y = pred),
            colour = "red", linewidth = 1.3) +
  
  # x 轴（0–60）
  scale_x_continuous(
    breaks = seq(0,60, by = 12),
    labels = seq(0,60, by = 12),
    name = "Time since surgery (months)"
  ) +
  
  # y 轴（效应量）
  scale_y_continuous(name = "Effect size (yi)") +
  
  # 坐标范围
  coord_cartesian(xlim = c(0, 60)) +
  
  # 气泡大小
  scale_size(range = c(1, 6), name = "Participants (n)") +
  
  # 主题
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "plain"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_line(linewidth = rel(0.5), linetype = 2)
  ) +
  
  # 👇 注释
  annotate("text", x = 55, y = -1.2,
           label = paste0("K = ", total_k, " (", total_studies, " studies)"),
           hjust = 1, size = 3.2) +
  annotate("text", x = 50, y = -1.6,
           label = paste0("N = ", total_n),
           hjust = 1, size = 3.2)

#可视化3-knot RCS 

# 构造新的预测时间点
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 用 rcspline.eval 生成 3-knot RCS 基函数
new_rcs3 <- rcspline.eval(new_data$time_point, nk = 3, inclx = TRUE)

# 执行预测
pred <- predict(model_rcs3, newmods = new_rcs3)

# 整理结果
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 4️⃣ 计算总样本数、研究数
total_k <- model_rcs3$k
total_studies <- dat %>% filter(!is.na(vi)) %>% summarise(n = n_distinct(study_id)) %>% pull(n)
total_n <- dat %>%
  filter(!is.na(vi)) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = sum(n_total, na.rm = TRUE)) %>%
  pull(n)

# 5️⃣ 绘图
ggplot() +
  # 基准线
  geom_hline(yintercept = 0, colour = "dark grey") +
  
  # 每个 study_id 的线条
  geom_line(data = dat, aes(x = time_point, y = yi, group = study_id),
            alpha = 0.8, colour = "grey") +
  
  # 气泡散点
  geom_point(data = dat, aes(x = time_point, y = yi, size = n_post_exp),
             alpha = 0.3, colour = "black", fill = "grey") +
  
  # 外层更宽阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 3,
                                   ymax = pred + 3),
              fill = "black", alpha = 0.1) +
  
  # 内层窄阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 1,
                                   ymax = pred + 1),
              fill = "red", alpha = 0.2) +
  
  # 主曲线
  geom_line(data = new_data, aes(x = time_point, y = pred),
            colour = "red", linewidth = 1.3) +
  
  # x 轴（0–60）
  scale_x_continuous(
    breaks = seq(0,72, by = 12),
    labels = seq(0,72, by = 12),
    name = "Time since surgery (months)"
  ) +
  
  # y 轴（效应量）
  scale_y_continuous(name = "Effect size (yi)", limits = c(-10, 25)) +
  
  # 坐标范围
  coord_cartesian(xlim = c(0, 28)) +
  
  # 气泡大小
  scale_size(range = c(1, 6), name = "Participants (n)") +
  
  # 主题
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "plain"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_line(linewidth = rel(0.5), linetype = 2)
  ) +
  
  # 👇 注释
  annotate("text", x = 60, y = 2,
           label = paste0("K = ", total_k, " (", total_studies, " studies)"),
           hjust = 1, size = 3.2) +
  annotate("text", x = 55, y = 0.8,
           label = paste0("N = ", total_n),
           hjust = 1, size = 3.2)

#可视化4-Knot RCS Model
# 构造新时间点（用于预测）
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 构造 4-knot RCS 的样条基函数矩阵
new_rcs4 <- rcspline.eval(new_data$time_point, nk = 4, inclx = TRUE)

# 使用模型进行预测
pred <- predict(model_rcs4, newmods = new_rcs4)

# 整理预测结果到 new_data 中
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 4️⃣ 计算总样本数、研究数
total_k <- model_rcs4$k
total_studies <- dat %>% filter(!is.na(vi)) %>% summarise(n = n_distinct(study_id)) %>% pull(n)
total_n <- dat %>%
  filter(!is.na(vi)) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(n = sum(n_total, na.rm = TRUE)) %>%
  pull(n)

# 5️⃣ 绘图
ggplot() +
  # 基准线
  geom_hline(yintercept = 0, colour = "dark grey") +
  
  # 每个 study_id 的线条
  geom_line(data = dat, aes(x = time_point, y = yi, group = study_id),
            alpha = 0.8, colour = "grey") +
  
  # 气泡散点
  geom_point(data = dat, aes(x = time_point, y = yi, size = n_post_exp),
             alpha = 0.3, colour = "black", fill = "grey") +
  
  # 外层更宽阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 15,
                                   ymax = pred + 15),
              fill = "black", alpha = 0.1) +
  
  # 内层窄阴影
  geom_ribbon(data = new_data, aes(x = time_point,
                                   ymin = pred - 5,
                                   ymax = pred + 5),
              fill = "red", alpha = 0.2) +
  
  # 主曲线
  geom_line(data = new_data, aes(x = time_point, y = pred),
            colour = "red", linewidth = 1.3) +
  
  # x 轴（0–60）
  scale_x_continuous(
    breaks = seq(0,96, by = 12),
    labels = seq(0,96, by = 12),
    name = "Time since surgery (months)"
  ) +
  
  # y 轴（效应量）
  scale_y_continuous(name = "Effect size (yi)", limits = c(-50, 150)) +
  
  # 坐标范围
  coord_cartesian(xlim = c(0, 36)) +
  
  # 气泡大小
  scale_size(range = c(1, 6), name = "Participants (n)") +
  
  # 主题
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "plain"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_line(linewidth = rel(0.5), linetype = 2)
  ) +
  
  # 👇 注释
  annotate("text", x = 32, y = -32,
           label = paste0("K = ", total_k, " (", total_studies, " studies)"),
           hjust = 1, size = 3.2) +
  annotate("text", x = 30, y = -37,
           label = paste0("N = ", total_n),
           hjust = 1, size = 3.2)




#关键时间点预测
summarize_model_predictions <- function(model, time_seq = 1:120, time_labels = c(3, 6, 12,18,24,36)) {
  has_mods <- !is.null(model$call$mods)
  
  if (!has_mods) {
    preds <- predict(model)
    time_values <- time_seq
  } else {
    # 更稳妥的结构判断逻辑
    mods_formula <- deparse(model$call$mods)
    
    if (grepl("log\\(time_point\\)", mods_formula)) {
      preds <- predict(model, newmods = log(time_seq))
      time_values <- time_seq
    } else if (any(grepl("ns3_", names(model$data)))) {
      new_ns3 <- ns(time_seq, df = 3)
      colnames(new_ns3) <- paste0("ns3_", 1:ncol(new_ns3))
      preds <- predict(model, newmods = new_ns3)
      time_values <- time_seq
    } else if (any(grepl("ns4_", names(model$data)))) {
      new_ns4 <- ns(time_seq, df = 4)
      colnames(new_ns4) <- paste0("ns4_", 1:ncol(new_ns4))
      preds <- predict(model, newmods = new_ns4)
      time_values <- time_seq
    } else if (grepl("time_point", mods_formula)) {
      preds <- predict(model, newmods = time_seq)
      time_values <- time_seq
    } else {
      stop("Unsupported model structure")
    }
  }
  
  # 创建预测结果数据框
  pred_df <- data.frame(time = time_values, pred = preds$pred, ci.lb = preds$ci.lb, ci.ub = preds$ci.ub)
  
  # 关键时间点预测值
  key_points <- pred_df %>%
    filter(round(time) %in% time_labels) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(label = paste0(time, " mo: ", pred, " [", ci.lb, ", ", ci.ub, "]")) %>%
    select(time, label)
  
  # zero-crossing
  zero_cross <- NA
  for (i in seq_len(nrow(pred_df) - 1)) {
    if (pred_df$ci.ub[i] < 0 & pred_df$ci.ub[i+1] >= 0) {
      zero_cross <- approx(x = pred_df$ci.ub[i:(i+1)], y = pred_df$time[i:(i+1)], xout = 0)$y
      break
    }
  }
  
  # 最后时间点
  last_data <- max(model$data$time_point, na.rm = TRUE)
  
  list(
    key_predictions = key_points,
    zero_crossing = ifelse(is.na(zero_cross), "Not crossed", paste0(round(zero_cross, 1), " mo")),
    last_data_point = paste0(round(last_data, 1), " mo")
  )
}
summarize_model_predictions(model_linear)   # 对 linear 模型
summarize_model_predictions(model_log)      # 对 log 模型
summarize_model_predictions(model_ns3)      # 对样条-3 模型
summarize_model_predictions(model_ns4)      # 对样条-4 模型


# 预测图可视化
summarize_model_predictions <- function(model, time_seq = 1:30, time_labels = c(3,6,12,18,24)) {
  # 提取 mods 模型表达式作为字符
  mod_str <- paste(as.character(model$call$mods), collapse = " ")
  
  # 区分模型结构并生成预测
  if (is.null(model$call$mods)) {
    preds <- predict(model)  # 空模型（无协变量）
    time_values <- time_seq
  } else if (grepl("log", mod_str)) {
    preds <- predict(model, newmods = log(time_seq))
    time_values <- time_seq
  } else if (grepl("time_point", mod_str)) {
    preds <- predict(model, newmods = time_seq)
    time_values <- time_seq
  } else if (grepl("ns3_", mod_str)) {
    new_ns3 <- ns(time_seq, df = 3)
    colnames(new_ns3) <- paste0("ns3_", 1:ncol(new_ns3))
    preds <- predict(model, newmods = new_ns3)
    time_values <- time_seq
  } else if (grepl("ns4_", mod_str)) {
    new_ns4 <- ns(time_seq, df = 4)
    colnames(new_ns4) <- paste0("ns4_", 1:ncol(new_ns4))
    preds <- predict(model, newmods = new_ns4)
    time_values <- time_seq
  } else {
    stop("Unsupported model structure")
  }
  
  # 整理预测数据
  pred_df <- data.frame(time = time_values, pred = preds$pred, ci.lb = preds$ci.lb, ci.ub = preds$ci.ub)
  
  key_points <- pred_df %>%
    filter(round(time) %in% time_labels) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    mutate(label = paste0(time, " month ")) %>%
    select(time, pred, label)
  
  # zero-crossing（上限穿过0）
  zero_cross <- NA
  for (i in seq_len(nrow(pred_df) - 1)) {
    if (pred_df$ci.ub[i] < 0 & pred_df$ci.ub[i + 1] >= 0) {
      zero_cross <- approx(x = pred_df$ci.ub[i:(i + 1)], y = pred_df$time[i:(i + 1)], xout = 0)$y
      break
    }
  }
  
  last_data <- max(model$data$time_point, na.rm = TRUE)
  
  # 绘图
  p <- ggplot(pred_df, aes(x = time, y = pred)) +
    # 🔴 新增“跟随曲线的阴影层” (可自行修改带宽，比如 ±0.1)
    # 🔴 第一层阴影（较宽）
    geom_ribbon(aes(ymin = pred - 0.35, ymax = pred + 0.35),
                fill = "#B2C7E9", alpha = 0.4) +
    geom_ribbon(aes(ymin = pred - 0.1, ymax = pred + 0.1),
                fill = "#F9A1B2", alpha = 0.3) +
    geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "white", alpha = 0.25) +
    geom_line(color = "#D7263D", linewidth = 1.2) +
    geom_point(data = key_points, aes(x = time, y = pred), color = "grey60", size = 3, inherit.aes = FALSE) +
    geom_text(data = key_points, aes(x = time, y = pred, label = label),
              inherit.aes = FALSE, hjust = -0.2, vjust = -0.2, size = 3.8) +
    labs(
      title = "",
      x = "Time (months)",
      y = "Predicted effect size"
    ) +
    theme_minimal()+
    scale_y_continuous(
      name = "Predicted effect size",
      limits = c(-2, 0.5)   # 直接限制范围
    )+
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),  # 主网格线淡化
    panel.grid.minor = element_blank()                               # 取消次网格线
  )
  
  print(p)
  
  return(list(
    key_predictions = key_points,
    zero_crossing = ifelse(is.na(zero_cross), "Not crossed", paste0(round(zero_cross, 1), " mo")),
    last_data_point = paste0(round(last_data, 1), " mo"),
    plot = p
  ))
}
result <- summarize_model_predictions(model_log)
result$plot            # 再次查看图（可在 R Markdown 中使用）
result$key_predictions # 关键时间点预测
result$zero_crossing   # 穿过0的时间


#（待选择）
#四个模型形对比图（一张图）
# 1. 为 dat 添加自然样条列（df = 3）
ns3 <- ns(dat$time_point, df = 3)
colnames(ns3) <- paste0("ns3_", 1:ncol(ns3))
dat_ns3 <- cbind(dat, ns3)

# 2. 同样添加 df = 4 的样条列
ns4 <- ns(dat$time_point, df = 4)
colnames(ns4) <- paste0("ns4_", 1:ncol(ns4))
dat_ns4 <- cbind(dat, ns4)

# 3. 构建四种模型
model_linear <- rma.mv(yi, vi, random = ~ time_point | study_id, data = dat)
model_log    <- rma.mv(yi, vi, random = ~ time_point | study_id, mods = ~ log(time_point), data = dat)
model_ns3    <- rma.mv(yi, vi, random = ~ time_point | study_id, mods = ~ ns3_1 + ns3_2 + ns3_3, data = dat_ns3)
model_ns4    <- rma.mv(yi, vi, random = ~ time_point | study_id, mods = ~ ns4_1 + ns4_2 + ns4_3 + ns4_4, data = dat_ns4)

# 4. 新时间点（用于预测）
time_seq <- seq(1, 90, length.out = 200)

# 构建新数据框（含 spline 基函数矩阵）
new_ns3 <- ns(time_seq, df = 3)
colnames(new_ns3) <- paste0("ns3_", 1:ncol(new_ns3))

new_ns4 <- ns(time_seq, df = 4)
colnames(new_ns4) <- paste0("ns4_", 1:ncol(new_ns4))

# 5. 执行预测
pred_linear <- predict(model_linear)  # ✅ 不带 newmods，因为是空模型
pred_log    <- predict(model_log,    newmods = log(time_seq))
pred_ns3    <- predict(model_ns3,    newmods = new_ns3)
pred_ns4    <- predict(model_ns4,    newmods = new_ns4)


# 6. 合并结果
model_plot_data <- rbind(
  data.frame(time = time_seq, pred = pred_linear$pred, ci.lb = pred_linear$ci.lb, ci.ub = pred_linear$ci.ub, model = "Linear"),
  data.frame(time = time_seq, pred = pred_log$pred,    ci.lb = pred_log$ci.lb,    ci.ub = pred_log$ci.ub,    model = "Log"),
  data.frame(time = time_seq, pred = pred_ns3$pred,    ci.lb = pred_ns3$ci.lb,    ci.ub = pred_ns3$ci.ub,    model = "Spline-3"),
  data.frame(time = time_seq, pred = pred_ns4$pred,    ci.lb = pred_ns4$ci.lb,    ci.ub = pred_ns4$ci.ub,    model = "Spline-4")
)

# 7. 可视化
ggplot(model_plot_data, aes(x = time, y = pred, color = model, fill = model)) +
  geom_point(data = dat, aes(x = time_point, y = yi), inherit.aes = FALSE, color = "black") +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.15, color = NA) +
  scale_color_manual(values = c(
    "Linear" = "#1b9e77",
    "Log" = "#d95f02",
    "Spline-3" = "#7570b3",
    "Spline-4" = "#e7298a"
  )) +
  scale_fill_manual(values = c(
    "Linear" = "#1b9e77",
    "Log" = "#d95f02",
    "Spline-3" = "#7570b3",
    "Spline-4" = "#e7298a"
  )) +
  labs(
    x = "Time (months)",
    y = "Effect size (yi)",
    title = "Comparison of Model-Predicted Trajectories"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())


######
ggplot() +
  geom_point(data = dat, aes(x = time_point, y = yi), color = "black") +
  geom_line(data = new_data, aes(x = time_point, y = pred), color = "#AA3A49") +
  geom_ribbon(data = new_data, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "#D0908F", alpha = 0.2) +
  labs(x = "Time (months)", y = "Effect size (yi)") +
  theme_minimal()

new_data_plot <- filter(new_data, time_point <= 36)

ggplot() +
  # 气泡图
  geom_point(data = filter(dat, time_point <= 36), 
             aes(x = time_point, y = yi, size = n_post_exp, color = yi), alpha = 0.6) +
  geom_line(data = new_data_plot, aes(x = time_point, y = pred), color = "#AA3A49", size = 1.2) +
  geom_ribbon(data = new_data_plot, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.2) +
  scale_size_continuous(range = c(1, 4), name = "Sample Size") +
  scale_color_gradient2(low = "#4575b4", mid = "gray60", high = "gray60", midpoint = 0,
                        name = "Effect Size") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), name = "Time (months)") +
  labs(y = "Effect size (yi)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", size = 0.3)

# 可视化 linear 模型
# 新时间点数据
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 进行逐点预测
pred <- predict(model_linear, newmods = new_data$time_point)

# 将预测结果赋值到 new_data
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 绘图
ggplot() +
  geom_point(data = dat, aes(x = time_point, y = yi), color = "black", size = 1.8) +
  geom_line(data = new_data, aes(x = time_point, y = pred), color = "#AA3A49", linewidth = 1.2) +
  geom_ribbon(data = new_data, aes(x = time_point, ymin = ci.lb, ymax = ci.ub),
              fill = "#D0908F", alpha = 0.25) +
  labs(
    title = "Linear Model Prediction with 95% CI",
    x = "Time (months)",
    y = "Effect size (yi)"
  ) +
  theme_minimal()

new_data_plot <- filter(new_data, time_point <= 36)

ggplot() +
  # 气泡图
  geom_point(data = filter(dat, time_point <= 36), 
             aes(x = time_point, y = yi, size = n_post_exp, color = yi), alpha = 0.6) +
  geom_line(data = new_data_plot, aes(x = time_point, y = pred), color = "#AA3A49", size = 1.2) +
  geom_ribbon(data = new_data_plot, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.2) +
  scale_size_continuous(range = c(1, 4), name = "Sample Size") +
  scale_color_gradient2(low = "#4575b4", mid = "gray60", high = "gray60", midpoint = 0,
                        name = "Effect Size") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), name = "Time (months)") +
  labs(y = "Effect size (yi)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", size = 0.3)


#可视化3-knot RCS 

# 构造新的预测时间点
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 用 rcspline.eval 生成 3-knot RCS 基函数
new_rcs3 <- rcspline.eval(new_data$time_point, nk = 3, inclx = TRUE)

# 执行预测
pred <- predict(model_spline, newmods = new_rcs3)

# 整理结果
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 可视化
ggplot() +
  geom_point(data = dat, aes(x = time_point, y = yi), color = "black") +
  geom_line(data = new_data, aes(x = time_point, y = pred), color = "#AA3A49") +
  geom_ribbon(data = new_data, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "#D0908F", alpha = 0.2) +
  labs(x = "Time (months)", y = "Effect size (yi)", 
       title = "3-knot RCS Model Prediction with 95% CI") +
  theme_minimal()

new_data_plot <- filter(new_data, time_point <= 36)

ggplot() +
  # 气泡图
  geom_point(data = filter(dat, time_point <= 36), 
             aes(x = time_point, y = yi, size = n_post_exp, color = yi), alpha = 0.6) +
  geom_line(data = new_data_plot, aes(x = time_point, y = pred), color = "#AA3A49", size = 1.2) +
  geom_ribbon(data = new_data_plot, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.2) +
  scale_size_continuous(range = c(1, 4), name = "Sample Size") +
  scale_color_gradient2(low = "#4575b4", mid = "gray60", high = "gray60", midpoint = 0,
                        name = "Effect Size") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), name = "Time (months)") +
  labs(y = "Effect size (yi)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", size = 0.3)


#可视化4-Knot RCS Model

# 构造新时间点（用于预测）
new_data <- data.frame(time_point = seq(1, 90, length.out = 100))

# 构造 4-knot RCS 的样条基函数矩阵
new_rcs4 <- rcspline.eval(new_data$time_point, nk = 4, inclx = TRUE)

# 使用模型进行预测
pred <- predict(model_rcs4, newmods = new_rcs4)

# 整理预测结果到 new_data 中
new_data$pred <- pred$pred
new_data$ci.lb <- pred$ci.lb
new_data$ci.ub <- pred$ci.ub

# 绘制图像
ggplot() +
  geom_point(data = dat, aes(x = time_point, y = yi), color = "black", size = 1.8) +
  geom_line(data = new_data, aes(x = time_point, y = pred), color = "#AA3A49", linewidth = 1.2) +
  geom_ribbon(data = new_data, aes(x = time_point, ymin = ci.lb, ymax = ci.ub),
              fill = "#D0908F", alpha = 0.2) +
  labs(
    title = "4-Knot RCS Model Prediction with 95% CI",
    x = "Time (months)",
    y = "Effect size (yi)"
  ) +
  theme_minimal()

new_data_plot <- filter(new_data, time_point <= 36)

ggplot() +
  # 气泡图
  geom_point(data = filter(dat, time_point <= 36), 
             aes(x = time_point, y = yi, size = n_post_exp, color = yi), alpha = 0.6) +
  geom_line(data = new_data_plot, aes(x = time_point, y = pred), color = "#AA3A49", size = 1.2) +
  geom_ribbon(data = new_data_plot, aes(x = time_point, ymin = ci.lb, ymax = ci.ub), 
              fill = "gray", alpha = 0.2) +
  scale_size_continuous(range = c(1, 4), name = "Sample Size") +
  scale_color_gradient2(low = "#4575b4", mid = "gray60", high = "gray60", midpoint = 0,
                        name = "Effect Size") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), name = "Time (months)") +
  labs(y = "Effect size (yi)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed", size = 0.3)
