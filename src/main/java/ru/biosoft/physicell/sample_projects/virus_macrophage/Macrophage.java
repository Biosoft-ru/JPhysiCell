package ru.biosoft.physicell.sample_projects.virus_macrophage;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;

public class Macrophage extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // bookkeeping 
        Microenvironment microenvironment = pCell.getMicroenvironment();
        int nVirus = microenvironment.findDensityIndex( "virus" );

        CellDefinition pMacrophage = CellDefinition.getCellDefinition( "macrophage" );

        // digest virus particles inside me 
        double implicit_Euler_constant = ( 1.0 + dt * pCell.custom_data.get( "virus_digestion_rate" ) );
        phenotype.molecular.internalized_total_substrates[nVirus] /= implicit_Euler_constant;

        // check for contact with a cell
        Cell pTestCell = null;
        List<Cell> neighbors = get_possible_neighbors( pCell );

        //	for( int n=0; n < pCell.cells_in_my_container().size() ; n++ )
        for( int n = 0; n < neighbors.size(); n++ )
        {
            pTestCell = neighbors.get( n );
            // if it is not me and not a macrophage 
            if( pTestCell != pCell && pTestCell.type != pMacrophage.type )
            {
                // calculate distance to the cell 
                double[] displacement = VectorUtil.newDiff( pTestCell.position, pCell.position );
                double distance = VectorUtil.norm( displacement );
                double max_distance = pCell.phenotype.geometry.radius + pTestCell.phenotype.geometry.radius;
                max_distance *= 1.1;

                // if it is not a macrophage, test for viral load 
                // if high viral load, eat it. 

                if( pTestCell.phenotype.molecular.internalized_total_substrates[nVirus] > pCell.custom_data
                        .get( "min_virion_detection_threshold" ) && distance < max_distance )
                {
                    System.out.println( "\t\tnom nom nom" );
                    pCell.ingestCell( pTestCell );
                }
            }
        }
    }

    public static List<Cell> get_possible_neighbors(Cell pCell)
    {
        List<Cell> neighbors = new ArrayList<>();

        // First check the neighbors in my current voxel
        //  std::vector<Cell>::iterator neighbor;
        //  std::vector<Cell>::iterator end =
        //      pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].end();
        //  for( Cell neighbor = pCell.get_container().agent_grid[pCell.get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
        //  { neighbors.push_back( neighbor ); }

        for( Cell neighbor : pCell.get_container().agent_grid.get( pCell.get_current_mechanics_voxel_index() ) )
        {
            neighbors.add( neighbor );
        }

        for( int ind : pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()] )
        {
            for( Cell neighbor : pCell.get_container().agent_grid.get( ind ) )
            {
                neighbors.add( neighbor );
            }
        }
        return neighbors;
    }
    //  for (Cell neighborVoxel: pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()])
    //  {
    //      
    //  }
    //  std::vector<int>::iterator neighbor_voxel_index;
    //  std::vector<int>::iterator neighbor_voxel_index_end = 
    //      pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].end();

    //  for( neighbor_voxel_index = 
    //      pCell.get_container().underlying_mesh.moore_connected_voxel_indices[pCell.get_current_mechanics_voxel_index()].begin();
    //      neighbor_voxel_index != neighbor_voxel_index_end; 
    //      ++neighbor_voxel_index )
    //  {
    //      if(!is_neighbor_voxel(pCell, pCell.get_container().underlying_mesh.voxels[pCell.get_current_mechanics_voxel_index()].center, pCell.get_container().underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
    //          continue;
    //      end = pCell.get_container().agent_grid[neighbor_voxel_index].end();
    //      for(neighbor = pCell.get_container().agent_grid[neighbor_voxel_index].begin();neighbor != end; ++neighbor)
    //      { neighbors.push_back( neighbor ); }
    //  }

    //  return neighbors; 
    //}
}