#include <cmath>
#include <iostream>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <util/Scenario.hpp>

int main()
{
  PLLRootedTree species("((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);", false);
  PLLUnrootedTree gene("(A:0.2,B:0.2,(C:0.2,D:0.2):0.6);", false);
  GeneSpeciesMapping mapping;
  mapping.fillFromGeneLabels(gene.getLabels());

  for (auto model: {RecModel::UndatedDL, RecModel::UndatedDTL}) {
    RecModelInfo options;
    options.model = model;
    options.originationStrategy = OriginationStrategy::UNIFORM;
    options.madRooting = true;
    options.pruneSpeciesTree = false;
    ReconciliationEvaluation evaluator(species, gene, mapping, options, "");
    constexpr double rate = 1e-6;
    evaluator.setRates(model == RecModel::UndatedDL
        ? Parameters(rate, rate) : Parameters(rate, rate, rate));
    double ll = evaluator.evaluate();
    if (!std::isfinite(ll) || ll > 0.0) {
      std::cerr << "Invalid MAD-weighted reconciliation log likelihood: " << ll << '\n';
      return 1;
    }
    Scenario scenario;
    evaluator.inferMLScenario(scenario);
    if (scenario.getEventCount(ReconciliationEventType::EVENT_S) != 3 ||
        scenario.getEventCount(ReconciliationEventType::EVENT_D) != 0) {
      std::cerr << "Unexpected reconciliation for a matching four-taxon topology\n";
      return 1;
    }
  }
}
