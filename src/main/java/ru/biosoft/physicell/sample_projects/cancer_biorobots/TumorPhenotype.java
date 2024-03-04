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

        double damageRate = SignalBehavior.getSingleSignal( pCell, "custom:damage_rate" );
        double repairRate = SignalBehavior.getSingleSignal( pCell, "custom:repair_rate" );
        double drugDeathRate = SignalBehavior.getSingleSignal( pCell, "custom:drug_death_rate" );

        double drug = SignalBehavior.getSingleSignal( pCell, "therapeutic" );

        double maxDamage = 1.0 * damageRate / ( 1e-16 + repairRate );

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

        //        double temp = drug;

        // reuse temp as much as possible to reduce memory allocations etc.
        //        temp *= dt;
        //        temp *= damageRate;

        damage += ( drug * dt * damageRate ) / ( repairRate * dt + 1 );//temp; // d_prev + dt*chemo*damage_rate

        //        temp = repairRate;
        //        temp *= dt;
        //        temp += 1.0;
        damage /= ( repairRate * dt + 1 );//temp; // (d_prev + dt*chemo*damage_rate)/(1 + dt*repair_rate)

        // then, see if the cell undergoes death from the therapy
        double temp = dt * damage * drugDeathRate / maxDamage;
//        temp = dt;
//        temp *= damage;
//        temp *= drugDeathRate;
//        temp /= maxDamage; // dt*(damage/max_damage)*death_rate

        // make sure we write the damage (not current a behavior)
        pCell.state.damage = damage;
        if( PhysiCellUtilities.UniformRandom() <= temp )
        {
            // pCell.start_death( apoptosis_model_index );
            SignalBehavior.setSingleBehavior( pCell, "apoptosis", 9e99 );
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;
        }
    }
}