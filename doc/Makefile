FILE = documentation

# Pasta para FILEs temporários
BUILD_DIR := build

# Cria a pasta build se não existir
$(shell mkdir -p $(BUILD_DIR))

PDF := $(FILE).pdf

all:
	@echo "Compilando $(FILE).tex..."
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) -aux-directory=$(BUILD_DIR) -include-directory=. $(FILE).tex
	-@if grep -q "Package biblatex" $(BUILD_DIR)/$(FILE).log 2>/dev/null; then \
		echo "Rodando biber..."; \
		biber $(BUILD_DIR)/$(FILE); \
	fi
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(FILE).tex
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(FILE).tex
	@echo "Movendo PDF para pasta principal..."
	@cp $(BUILD_DIR)/$(PDF) ./

clean:
	rm -rf $(BUILD_DIR)/*.aux $(BUILD_DIR)/*.log $(BUILD_DIR)/*.bbl \
		$(BUILD_DIR)/*.blg $(BUILD_DIR)/*.toc $(BUILD_DIR)/*.nav \
		$(BUILD_DIR)/*.out $(BUILD_DIR)/*.snm $(BUILD_DIR)/*.synctex.gz
	rm -f $(PDF)