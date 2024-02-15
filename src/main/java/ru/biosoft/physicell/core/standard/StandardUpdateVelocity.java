package ru.biosoft.physicell.core.standard;

import java.util.Set;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellContainer;
import ru.biosoft.physicell.core.CellFunctions;
import ru.biosoft.physicell.core.Phenotype;

public class StandardUpdateVelocity implements CellFunctions.UpdateVelocity
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        if( pCell.functions.add_cell_basement_membrane_interactions != null )
        {
            pCell.functions.add_cell_basement_membrane_interactions.execute( pCell, phenotype, dt );
        }


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
        Set<Cell> neighbors = container.agent_grid.get( pCell.get_current_mechanics_voxel_index() );
        for( Cell neighbor : neighbors )
        {
            pCell.addPotentials( neighbor );
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
        int voxelIndex = pCell.get_current_mechanics_voxel_index();
        double[] center = container.underlying_mesh.voxels[voxelIndex].center;
        int[] neighborVoxels = container.underlying_mesh.moore_connected_voxel_indices[voxelIndex];
        for( int neighbor_voxel_index : neighborVoxels )
        {
            if( !Cell.isNeighborVoxel( pCell, center, container.underlying_mesh.voxels[neighbor_voxel_index].center,
                    neighbor_voxel_index ) )
                continue;
            for( Cell neighbor : container.agent_grid.get( neighbor_voxel_index ) )
            {
                pCell.addPotentials( neighbor );
            }
        }
        pCell.updateMotilityVector( dt );
        VectorUtil.sum( pCell.velocity, phenotype.motility.motilityVector );
    }
}