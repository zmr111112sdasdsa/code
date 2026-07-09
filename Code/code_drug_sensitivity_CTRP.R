
# 假设你已经加载了数据并且数据框命名为 'data'
exp$APEX1
# 1. 数据整理
data <- data.frame(
  CellLine = c("ACHN", "OSRC2", "VMRCRCZ", "RCC10RGB","CAKI2", "KMRC20", "KMRC3", 
               "TUHR14TKB", "769P", "CAKI1", "CAL54", "TUHR10TKB", "VMRCRCW", "TUHR4TKB", 
               "SNU1272", "A498", "786O", "KMRC1", "KMRC2", "BFTC909"),
  APEX1 = c(6.986775, 6.983507, 6.535420, 6.579509, 7.449372, 6.573176, 7.268257, 
            6.521352, 6.811664, 7.222756, 7.326906, 6.720073, 6.487379, 6.377776, 
            7.488089, 6.809247,  7.030049, 6.524730, 6.856923, 6.742720),
  
  erastin = c(13.737, 14.112, NA, 10.560, 11.808, 11.694, 13.518, 9.7769, 12.794, 
              14.620, 10.571, NA, 7.4045, 11.163, 13.233, 12.121, 14.373, 11.097, 
              11.567, 14.423),
  sorafenib = c(15.321, 12.289, 10.388, 11.285, 13.018, 13.312, 12.669, 13.2060, 
                13.272, 14.292, 12.564, 12.937, 14.4850, 13.753, 13.094, 16.334, 
                13.669, 13.771, 12.682, 12.229)
)

# 2. 按 APEX1 中位数切分成两组（高表达组和低表达组）
median_APEX1 <- median(data$APEX1, na.rm = TRUE)
data$APEX1_group <- ifelse(data$APEX1 > median_APEX1, "High", "Low")

# 3. 对比高表达组和低表达组的 AUC 值（使用 erastin 和 sorafenib）
# 我们假设数据没有缺失值，或者我们已经处理了缺失数据

# erastin 对比
t_test_erastin <- t.test(erastin ~ APEX1_group, data = data)
# sorafenib 对比
t_test_sorafenib <- t.test(sorafenib ~ APEX1_group, data = data)

# 4. 输出对比结果
cat("Erastin AUC comparison between high and low APEX1 expression groups:\n")
print(t_test_erastin)

cat("\nSorafenib AUC comparison between high and low APEX1 expression groups:\n")
print(t_test_sorafenib)

# 5. 结果可视化（可选）
# 绘制箱线图（boxplot）显示 AUC 值的分布
library(ggplot2)

# Erastin
ggplot(data, aes(x = APEX1_group, y = erastin, fill = APEX1_group)) +
  geom_boxplot() +
  labs(title = "Erastin AUC by APEX1 Expression Group", y = "AUC", x = "APEX1 Expression Group") +
  theme_minimal()

# Sorafenib
ggplot(data, aes(x = APEX1_group, y = sorafenib, fill = APEX1_group)) +
  geom_boxplot() +
  labs(title = "Sorafenib AUC by APEX1 Expression Group", y = "AUC", x = "APEX1 Expression Group") +
  theme_minimal()
