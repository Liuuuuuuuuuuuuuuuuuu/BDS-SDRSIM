# =====================[ 編譯設定 ]=====================
CC      = gcc
CFLAGS = -I. -O2 -Wall -Wextra -fopenmp -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

OBJS = main.o globals.o bch.o navbits.o channel.o \
        bdssim.o rinex.o orbits.o coord.o path.o
	
bds-sim: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) -lm

tests/prn_test: tests/prn_test.o globals.o bch.o navbits.o channel.o rinex.o orbits.o coord.o path.o bdssim.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

tests/test_beidou: tests/test_beidou.c
	$(CC) $(CFLAGS) $< -o $@

tests/test_prn: tests/test_prn.o globals.o bch.o navbits.o channel.o rinex.o orbits.o coord.o path.o bdssim.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

tests/test_nh_prn: tests/test_nh_prn.o globals.o bch.o navbits.o channel.o rinex.o orbits.o coord.o path.o bdssim.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

tests/test_orbits: tests/test_orbits.o globals.o bch.o navbits.o channel.o rinex.o orbits.o coord.o path.o bdssim.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

tests/test_timeconv: tests/test_timeconv.c
	$(CC) $(CFLAGS) $< -o $@

check: tests/test_prn tests/test_nh_prn
	./tests/test_prn
	./tests/test_nh_prn

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o bds-sim tests/prn_test tests/test_beidou \
        tests/test_prn tests/test_nh_prn tests/test_nh_prn.o \
        tests/test_orbits tests/test_orbits.o \
        tests/test_timeconv *.bin *.gz
	@echo "已清除中間檔與 .gz 暫存檔"

# =====================[ 下載區 ]=====================
# 使用範例：
#   make download-brdm20251760          (下載 2025/176 BRDM)
#   make download-brdc0010.24n          (下載 2024/001 BRDC)
#
# 命名規則：
#   BRDM ：brdmYYYYDDD0          → 自動補上完整 MGEX 檔名
#   BRDC ：brdcDDD0.YYn          → 與 IGS 老檔案一致
#
download-%:
	@NAME_LOWER=$*; \
	NAME=$$(echo $$NAME_LOWER | tr a-z A-Z); \
	if echo $$NAME | grep -q '^BRDM'; then \
		# ------------- BRDM (多GNSS) ------------- \
		YEAR=$$(echo $$NAME | cut -c5-8); \
		DOY=$$(echo $$NAME | cut -c9-11); \
		FILE=BRDM00DLR_S_$${YEAR}$${DOY}0000_01D_MN.rnx.gz; \
		SUBDIR=brdm; \
	elif echo $$NAME | grep -q '^BRDC'; then \
		# ------------- BRDC (GPS/GLONASS) -------- \
		YEAR=20$$(echo $$NAME | sed 's/.*\.\([0-9][0-9]\)N/\1/' ); \
		DOY=$$(echo $$NAME | cut -c5-7); \
		FILE=$${NAME}.gz; \
		SUBDIR=brdc; \
	else \
		echo "無法辨識檔名 $$NAME_LOWER"; exit 1; \
	fi; \
	URL=https://cddis.nasa.gov/archive/gnss/data/daily/$${YEAR}/brdc/$$FILE; \
	echo "下載 $$URL"; \
	curl -s -L -n -c cookies.txt -b cookies.txt \
	     -H "User-Agent: Mozilla/5.0" \
	     -o $$FILE $$URL; \
	if grep -q '<html' $$FILE; then \
		echo "下載失敗，收到 HTML 錯誤頁 (未授權?)"; \
		rm -f $$FILE; exit 1; \
	fi; \
	if file $$FILE | grep -q 'gzip compressed'; then \
		gunzip -f $$FILE; \
		echo "解壓完成 → $${FILE%.gz}"; \
	else \
		echo "檔案已是解壓狀態：$$FILE"; \
	fi

.PHONY: all clean download-%

