package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.standard.O2based;

public class TumorPhenotype extends UpdatePhenotype
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        double damage = SignalBehavior.getSingleSignal( pCell, "damage" );

        double damage_rate = SignalBehavior.getSingleSignal( pCell, "custom:damage_rate" );
        double repair_rate = SignalBehavior.getSingleSignal( pCell, "custom:repair_rate" );
        double drug_death_rate = SignalBehavior.getSingleSignal( pCell, "custom:drug_death_rate" );

        double drug = SignalBehavior.getSingleSignal( pCell, "therapeutic" );

        double max_damage = 1.0 * damage_rate / ( 1e-16 + repair_rate );

        // if I'm dead, don't bother. disable my phenotype rule
        if( SignalBehavior.getSingleSignal( pCell, "dead" ) > 0.5 )
        {
            pCell.functions.updatePhenotype = null;
            return;
        }

        // first, vary the cell birth and death rates with oxygenation

        // std::cout << get_single_behavior( pCell , "cycle entry") << " vs ";
        new O2based().execute( pCell, phenotype, dt );
        // std::cout << get_single_behavior( pCell , "cycle entry") << std::endl;

        // the update the cell damage

        // dD/dt = alpha*c - beta-D by implicit scheme

        double temp = drug;

        // reuse temp as much as possible to reduce memory allocations etc.
        temp *= dt;
        temp *= damage_rate;

        damage += temp; // d_prev + dt*chemo*damage_rate

        temp = repair_rate;
        temp *= dt;
        temp += 1.0;
        damage /= temp; // (d_prev + dt*chemo*damage_rate)/(1 + dt*repair_rate)

        // then, see if the cell undergoes death from the therapy
        temp = dt;
        temp *= damage;
        temp *= drug_death_rate;
        temp /= max_damage; // dt*(damage/max_damage)*death_rate

        // make sure we write the damage (not current a behavior)
        pCell.state.damage = damage;
        if( damage > 0 )
        {
            System.out.println( damage );
        }
        if( PhysiCellUtilities.UniformRandom() <= temp )
        {
            // pCell.start_death( apoptosis_model_index );
            SignalBehavior.setSingleBehavior( pCell, "apoptosis", 9e99 );
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;
        }
    }
}