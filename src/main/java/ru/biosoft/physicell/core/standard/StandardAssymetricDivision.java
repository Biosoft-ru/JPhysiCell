package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.AsymmetricDivision;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CellDivision;
import ru.biosoft.physicell.core.RandomGenerator;

public class StandardAssymetricDivision extends CellDivision
{
    RandomGenerator generator;
    public StandardAssymetricDivision(RandomGenerator generator)
    {
        this.generator = generator;
    }
    
    public void execute( Cell parent, Cell daughter ) throws Exception
    {
        
       // double total = pCell_parent->phenotype.cycle.asymmetric_division.probabilities_total();
        AsymmetricDivision parentDivision = parent.phenotype.cycle.getAsymmetricDivision();
        double total = parentDivision.getTotalProbability();
        if (total > 1.0)
        {
            double sym_div_prob = parentDivision.getProbability( parent.typeName ) + 1.0 - total;
//            double sym_div_prob = parentDefinition.cycle.asymmetric_division_probabilities[pCell_parent->type] + 1.0 - total;
            if (sym_div_prob < 0.0)
            { 
                throw new Exception("Error: Asymmetric division probabilities for " + parent.typeName + " sum to greater than 1.0 and cannot be normalized.");
//                throw //std::runtime_error("Error: Asymmetric division probabilities for " + pCD_parent->name + " sum to greater than 1.0 and cannot be normalized.");
            }
            parentDivision.setProbability(parent.typeName, sym_div_prob);
            daughter.phenotype.cycle.getAsymmetricDivision().setProbability(daughter.typeName, sym_div_prob);
        }
        double r = generator.UniformRandom();
        double[] probabilities = parentDivision.getProbabilities();
        for( int i=0; i < probabilities.length; i++ )
        {
            if( r <= probabilities[i] )
            {
                if (i != daughter.type) // only convert if the daughter is not already the correct type
                {
                    daughter.convert( daughter.getModel().getCellDefinition( i) ); 
                }
                return;
            } 
            r -= probabilities[i];
        }
        // if we're here, then do not do asym div
        return;
    }
    
//    void standard_asymmetric_division_function( Cell* pCell_parent, Cell* pCell_daughter )
//    {
//        Cell_Definition* pCD_parent = cell_definitions_by_name[pCell_parent->type_name];
//        double total = pCell_parent->phenotype.cycle.asymmetric_division.probabilities_total();
//        if (total > 1.0)
//        {
//            double sym_div_prob = pCell_parent->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[pCell_parent->type] + 1.0 - total;
//            if (sym_div_prob < 0.0)
//            { 
//                throw std::runtime_error("Error: Asymmetric division probabilities for " + pCD_parent->name + " sum to greater than 1.0 and cannot be normalized.");
//            }
//            pCell_parent->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[pCell_parent->type] = sym_div_prob;
//            pCell_daughter->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[pCell_daughter->type] = sym_div_prob;
//        }
//        double r = UniformRandom();
//        for( int i=0; i < pCD_parent->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.size(); i++ )
//        {
//            if( r <= pCell_parent->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[i] )
//            {
//                if (i != pCell_daughter->type) // only convert if the daughter is not already the correct type
//                { pCell_daughter->convert_to_cell_definition( *cell_definitions_by_index[i] ); }
//                return;
//            }
//            r -= pCell_parent->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities[i];
//        }
//        // if we're here, then do not do asym div
//        return;
//    }
}
