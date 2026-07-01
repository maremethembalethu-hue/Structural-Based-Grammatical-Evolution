JAVAC = javac
JAVA = java

MAIN_SBGE = EnergyPredictSBGE

SRC = ExpressionEvaluator.java EnergyPredictSBGE.java GEIndividual.java BNFGrammar.java

compile:
	$(JAVAC) $(SRC)

run-sbge:
	@echo Running Structural-Based GE (EnergyPredictSBGE)...
	$(JAVA) $(MAIN_SBGE)

clean:
	rm -f *.class