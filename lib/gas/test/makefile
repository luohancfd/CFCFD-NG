TEST_PROG := ltest.lua

RESULT_FILES := user-defined-gas-model-test.result \
	composite-gas-model-test.result \
	generalised-Arrhenius-test.result

test: prepare $(RESULT_FILES)
	./$(TEST_PROG) user-defined-gas-model-test.result
	./$(TEST_PROG) composite-gas-model-test.result
	./$(TEST_PROG) generalised-Arrhenius-test.result
prepare: 
	- rm $(RESULT_FILES)

clean: prepare
	- rm *.x
	- rm $(RESULT_FILES)

user-defined-gas-model-test.result :
	./user-defined-gas-model-test.x

composite-gas-model-test.result :
	./composite-gas-model-test.x

generalised-Arrhenius-test.result :
	./generalised-Arrhenius-test.x

