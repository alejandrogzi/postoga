# Makefile

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.9.3-devel

.PHONY: configure test create --conda

PYTHON_SCRIPT=postoga_test.py
PROJECT := $(CURDIR)/POSTOGA_TEST
CONFIG=configure.sh
VERSION="0.9.3-devel"
TEST_PATH := $(CURDIR)/supply/test
THRESHOLD := 0.95
CLASS := I,PI
FMT := gtf

configure:
	@echo ""
	@echo "> postoga v.$(VERSION) configuration!"
	@echo ">> this may take some minutes..."
	@echo ""
	$(eval CONDA_FLAG := $(filter --conda,$(MAKECMDGOALS)))
	bash $(CONFIG) $(CONDA_FLAG)

--conda:
	@:

create:
	@echo "#!/usr/bin/env python3" > $(PYTHON_SCRIPT)
	@echo "from run import parser, TogaDir" >> $(PYTHON_SCRIPT)
	@echo "" >> $(PYTHON_SCRIPT)
	@echo "def main():" >> $(PYTHON_SCRIPT)
	@echo "    args = parser()" >> $(PYTHON_SCRIPT)
	@echo "    TogaDir(args).run()" >> $(PYTHON_SCRIPT)
	@echo "" >> $(PYTHON_SCRIPT)
	@echo "if __name__ == \"__main__\":" >> $(PYTHON_SCRIPT)
	@echo "    main()" >> $(PYTHON_SCRIPT)


test: create
	@python3 $(PYTHON_SCRIPT) base --togadir $(TEST_PATH) --outdir $(PROJECT) -th $(THRESHOLD) -bc $(CLASS) -to $(FMT) --extract || (echo "> ERROR: test failed!" && exit 1)
	@rm $(PYTHON_SCRIPT)
	@mv $(PROJECT)'/postoga.log' './_postoga.test.log'
	@rm -r $(PROJECT)
	@echo ""
	@echo "> test completed!"
