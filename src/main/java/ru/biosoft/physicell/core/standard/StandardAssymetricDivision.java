package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.AsymmetricDivision;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CellDivision;

public class StandardAssymetricDivision extends CellDivision
{
    public void execute( Cell parent, Cell daughter ) throws Exception
    {       
        AsymmetricDivision parentDivision = parent.phenotype.cycle.getAsymmetricDivision();
        double total = parentDivision.getTotalProbability();
        if (total > 1.0)
        {
            double sym_div_prob = parentDivision.getProbability( parent.typeName ) + 1.0 - total;
            if (sym_div_prob < 0.0)
            { 
                throw new Exception("Error: Asymmetric division probabilities for " + parent.typeName + " sum to greater than 1.0 and cannot be normalized.");
            }
            parentDivision.setProbability(parent.typeName, sym_div_prob);
            daughter.phenotype.cycle.getAsymmetricDivision().setProbability(daughter.typeName, sym_div_prob);
        }
        double r = parent.getModel().getRNG().UniformRandom();
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
    }
    
    @Override
    public String getName()
    {
        return "Standard asymmetric division";
    }
}
