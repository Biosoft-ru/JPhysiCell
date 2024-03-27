package ru.biosoft.physicell.sample_projects.ecoli_acetic_switch;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.VolumeUpdate;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.fba.IntracellularFBA;

public class update_cell extends VolumeUpdate
{
    double O2_Km;
    double O2_Vmax;
    double glc_Km;
    double glc_Vmax;
    double lac_Km;
    double lac_Vmax;

    public update_cell(Model model)
    {
        O2_Km = model.getParameterDouble( "oxygen_Km" );
        O2_Vmax = model.getParameterDouble( "oxygen_Vmax" );
        glc_Km = model.getParameterDouble( "glucose_Km" );
        glc_Vmax = model.getParameterDouble( "glucose_Vmax" );
        lac_Km = model.getParameterDouble( "acetate_Km" );
        lac_Vmax = model.getParameterDouble( "acetate_Vmax" );
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        Microenvironment m = pCell.getMicroenvironment();
        IntracellularFBA intra = (IntracellularFBA)phenotype.intracellular;
        phenotype.secretion.setSecretionToZero();
        phenotype.secretion.setUptakeToZero();

        String oxygen_name = "oxygen";
        String glucose_name = "glucose";
        String acetate_name = "acetate";

        String oxygen_flux_id = intra.exchange_flux_density_map( oxygen_name );
        String glucose_flux_id = intra.exchange_flux_density_map( glucose_name );
        String acetate_flux_id = intra.exchange_flux_density_map( acetate_name );

        int oxygen_idx = m.findDensityIndex( oxygen_name );
        int glucose_idx = m.findDensityIndex( glucose_name );
        int acetate_idx = m.findDensityIndex( acetate_name );

        int voxelIndex = pCell.currentVoxelIndex;
        double[] density = m.nearestDensity( voxelIndex );
        double oxygen_density = density[oxygen_idx];
        double glucose_density = density[glucose_idx]; // dived by voxel size?
        double acetate_density = density[acetate_idx]; // dived by voxel size?

        double oxygen_flux_bound = -1 * ( O2_Vmax * oxygen_density ) / ( oxygen_density + O2_Km );
        intra.model.setReactionLowerBound( oxygen_flux_id, oxygen_flux_bound );
        double glucose_flux_bound = -1 * ( glc_Vmax * glucose_density ) / ( glucose_density + glc_Km );
        intra.model.setReactionLowerBound( glucose_flux_id, glucose_flux_bound );
        double acetate_flux_bound = -1 * ( lac_Vmax * acetate_density ) / ( acetate_density + lac_Km );
        intra.model.setReactionLowerBound( acetate_flux_id, acetate_flux_bound );

        System.out.println( "Oxygen density: " + oxygen_density + " => " + intra.model.fbcModel.getLowerBound( oxygen_flux_id ) );
        System.out.println( "Glucose density: " + glucose_density + " => " + intra.model.fbcModel.getLowerBound( glucose_flux_id ) );
        System.out.println( "Acetate density: " + acetate_density + " => " + intra.model.fbcModel.getLowerBound( acetate_flux_id ) );

        intra.model.runFBA();

        //                System.out.println( "Optimized R_EX_o2_e " + intra.model.getFlux( "R_EX_o2_e" ) );
        if( intra.model.getSolutionStatus() )
        {
            double oxygen_flux = intra.model.getFlux( oxygen_flux_id );
            //            System.out.println( "Oxygen: " + oxygen_flux );
            //	std::cout << "Oxygen flux: " << oxygen_flux << std::endl;
            double glucose_flux = intra.model.getFlux( glucose_flux_id );
            //            System.out.println( "Glucose: " + glucose_flux );
            //	std::cout << "glucose flux: " << glucose_flux << std::endl;
            double acetate_flux = intra.model.getFlux( acetate_flux_id );
            //            System.out.println( "Acetate: " + acetate_flux );
            //	std::cout << "acetate flux: " << acetate_flux << std::endl;

            if( oxygen_flux < 0 )
                phenotype.secretion.uptakeRates[oxygen_idx] = Math.abs( oxygen_flux / oxygen_density );

            if( glucose_flux < 0 )
                phenotype.secretion.uptakeRates[glucose_idx] = Math.abs( glucose_flux / glucose_density );

            if( acetate_flux < 0 )
                phenotype.secretion.uptakeRates[acetate_idx] = Math.abs( acetate_flux / acetate_density );

            else if( acetate_flux > 0 )
                phenotype.secretion.secretionRates[acetate_idx] = Math.abs( acetate_flux / acetate_density );

        }
        else
        {
            System.out.println( "Apoptosis" );
            int nApoptosis = phenotype.death.findDeathModelIndex( "apoptosis" );
            pCell.startDeath( nApoptosis );
            //            System.out.println( "Could not find feasible solution" );

            // Check energy production to see is the cell is able
            // to reach a threshold. Otherwise enter apoptosis

        }
    }
}