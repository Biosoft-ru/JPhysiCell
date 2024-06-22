package ru.biosoft.physicell.core.standard;

import java.util.Set;

import ru.biosoft.physicell.biofvm.CartesianMesh;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdateVelocity;
import ru.biosoft.physicell.core.Mechanics;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;

public class StandardUpdateVelocity extends UpdateVelocity
{
    double[][] repulsion;
    double[][] adhesion;

    double simple_pressure_scale = 0.027288820670331;

    public static double addPTime = 0;
    public static double totalTime = 0;
    boolean isInit = false;

    public void init(Model model)
    {
        int defs = model.getDefinitionsCount();
        repulsion = new double[defs][defs];
        adhesion = new double[defs][defs];
        for( int i = 0; i < defs; i++ )
        {
            for( int j = 0; j < defs; j++ )
            {
                CellDefinition cdi = model.getCellDefinition( i );
                CellDefinition cdj = model.getCellDefinition( j );

                Mechanics mi = cdi.phenotype.mechanics;
                Mechanics mj = cdj.phenotype.mechanics;

                repulsion[i][j] = Math.sqrt( mi.cellCellRepulsionStrength * mj.cellCellRepulsionStrength );

                double adhesion_ii = mi.cellCellAdhesionStrength * mj.cellAdhesionAffinities[j];
                double adhesion_jj = mj.cellCellAdhesionStrength * mi.cellAdhesionAffinities[i];
                adhesion[i][j] = Math.sqrt( adhesion_ii * adhesion_jj );

            }
        }
        isInit = true;
    }
    
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
//        double start = System.nanoTime();
        if( pCell.functions.membraneInteraction != null )
        {
            pCell.functions.membraneInteraction.execute( pCell, phenotype, dt );
        }
        if( !isInit )
            init( pCell.getModel() );

        pCell.state.simplePressure = 0.0;
        pCell.state.neighbors.clear(); // new 1.8.0

        //First check the neighbors in my current voxel
        //        std::vector<Cell*>::iterator neighbor;
        //        std::vector<Cell*>::iterator end = pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].end();
        //        for(neighbor = pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
        //        {
        //            pCell.add_potentials(*neighbor);
        //        }
        CellContainer container = pCell.get_container();
        Set<Cell> neighbors = container.agentGrid.get( pCell.get_current_mechanics_voxel_index() );
        for( Cell neighbor : neighbors )
        {
            //                     
//            double start2 = System.nanoTime();
//            pCell.addPotentials( neighbor );
            addPotentials2( pCell, neighbor );
//            addPTime += System.nanoTime() - start2;
            
            //            if (p1 != p2)
            //            System.out.println( "P1 "+p1+" | P2 "+p2 );
        }
        //        std::vector<int>::iterator neighbor_voxel_index;
        //        std::vector<int>::iterator neighbor_voxel_index_end = 
        //            pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].end();               
        //        for( neighbor_voxel_index = 
        //            pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].begin();
        //            neighbor_voxel_index != neighbor_voxel_index_end; 
        //            ++neighbor_voxel_index )
        //        {
        //            if(!is_neighbor_voxel(pCell, pCell.get_container().underlying_mesh.voxels[pCell.get_current_mechanics_voxel_index()].center, pCell.get_container().underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
        //                continue;
        //            end = pCell.get_container().agent_grid[*neighbor_voxel_index].end();
        //            for(neighbor = pCell.get_container().agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
        //            {
        //                pCell.add_potentials(neighbor);
        //            }
        //        }
        CartesianMesh mesh = container.mesh;
        int voxelIndex = pCell.get_current_mechanics_voxel_index();
        double[] center = mesh.voxels[voxelIndex].center;
        int[] neighborVoxels = mesh.moore_connected_voxel_indices[voxelIndex];
        for( int neighborIndex : neighborVoxels )
        {
            if( !Cell.isNeighborVoxel( pCell, center, mesh.voxels[neighborIndex].center, neighborIndex ) )
            {
                continue;
            }
            for( Cell neighbor : container.agentGrid.get( neighborIndex ) )
            {
                //                                
//                double start2 = System.nanoTime();
//                pCell.addPotentials( neighbor );
                addPotentials2( pCell, neighbor );
//                addPTime += System.nanoTime() - start2;
                //                if (p1 != p2)
                //                System.out.println( "P1 "+p1+" | P2 "+p2 );
            }
        }
        pCell.updateMotilityVector( dt );
        VectorUtil.sum( pCell.velocity, phenotype.motility.motilityVector );
//        totalTime += System.nanoTime() - start;
    }

    public double addPotentials2(Cell c1, Cell c2)
    {
        if( c1.ID == c2.ID )
            return 0;
            
        int i = c1.type;
        int j = c2.type;
        double r1 = c1.phenotype.geometry.radius;
        double r2 = c2.phenotype.geometry.radius;
        double distance = 0;
        double[] displacement = new double[3];
        for( int k = 0; k < 3; k++ )
        {
            displacement[k] = c1.position[k] - c2.position[k];
            distance += displacement[k] * displacement[k];
        }

        distance = Math.max( Math.sqrt( distance ), 0.00001 );
        double R = r1 + r2;
        double temp_r;
        if( distance > R )
            temp_r = 0;
        else
        {
            temp_r = 1 - distance / R;
            temp_r *= temp_r;
            c1.state.simplePressure += ( temp_r / simple_pressure_scale );
            temp_r *= repulsion[i][j];
        }

        double max_interactive_distance = c1.phenotype.mechanics.relMaxAdhesionDistance * r1
                + c2.phenotype.mechanics.relMaxAdhesionDistance * r2;

        if( distance < max_interactive_distance )
        {
            double temp_a = 1 - distance / max_interactive_distance;
            temp_a *= temp_a;
            temp_a *= adhesion[i][j];
            temp_r -= temp_a;
            c1.state.neighbors.add( c2 );
        }
        if( Math.abs( temp_r ) < 1e-16 )
            return 0;
        temp_r /= distance;
        VectorUtil.axpy( c1.velocity, temp_r, displacement );
        return temp_r;
    }

    @Override
    public String display()
    {
        return "Standard velocity: cell-cell adhesion + biased motility";
    }

    @Override
    public String getName()
    {
        return "Standard velocity";
    }
}