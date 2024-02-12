package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.core.StandardModels;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;

public class TumorPhenotype extends update_phenotype
{
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        double damage = SignalBehavior.get_single_signal( pCell, "damage" );

        double damage_rate = SignalBehavior.get_single_signal( pCell, "custom:damage_rate" );
        double repair_rate = SignalBehavior.get_single_signal( pCell, "custom:repair_rate" );
        double drug_death_rate = SignalBehavior.get_single_signal( pCell, "custom:drug_death_rate" );

        double drug = SignalBehavior.get_single_signal( pCell, "therapeutic" );

        double max_damage = 1.0 * damage_rate / ( 1e-16 + repair_rate );

        // if I'm dead, don't bother. disable my phenotype rule
        if( SignalBehavior.get_single_signal( pCell, "dead" ) > 0.5 )
        {
            pCell.functions.updatePhenotype = null;
            return;
        }

        // first, vary the cell birth and death rates with oxygenation

        // std::cout << get_single_behavior( pCell , "cycle entry") << " vs ";
        new StandardModels.update_cell_and_death_parameters_O2_based().execute( pCell, phenotype, dt );
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
            pCell.functions.custom_cell_rule = null;
        }
    }
}