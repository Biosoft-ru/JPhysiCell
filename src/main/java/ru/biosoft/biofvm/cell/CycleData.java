package ru.biosoft.biofvm.cell;

import java.util.ArrayList;
import java.util.List;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
public class CycleData
{
    // this maps the end_phase_index to the link index in each 
    // phase_links[i]
    // So, index_inverse_map[i][j] = k, corresponds to 
    // phases[i], phase_links[i][k] (which links from phase i to phase j)
    // transition_rates[i][k] (the transition rate from phase i to phase j) 
    //        std::vector< std::unordered_map<int,int> > inverse_index_maps; 
    //    List<Map<Integer, Integer>> inverse_index_maps;

    private CycleModel pCycle_Model;

    String time_units;

    //        std::vector< std::vector<double> > transition_rates; 
    List<List<Double>> transition_rates;

    public int current_phase_index;
    public double elapsed_time_in_phase;

    public CycleData(CycleModel model)
    {
        //        inverse_index_maps = new ArrayList<>();//.resize(0); 
        pCycle_Model = model;
        time_units = "min";
        transition_rates = new ArrayList<>();//.resize( 0 );
        current_phase_index = 0;
        elapsed_time_in_phase = 0.0;
    }

    //    private void resizeListMap(List<Map<Integer, Integer>> list, int n)
    //    {
    //        if( n < list.size() )
    //        {
    //            for( int i = n; i < list.size(); i++ )
    //                list.remove( i );
    //        }
    //        else if( n > list.size() )
    //        {
    //            for( int i = 0; i < n - list.size(); i++ )
    //                list.add( new HashMap<Integer, Integer>() );
    //        }
    //    }

    private void resizeListList(List<List<Double>> list, int n)
    {
        if( n < list.size() )
        {
            for( int i = n; i < list.size(); i++ )
                list.remove( i );
        }
        else if( n > list.size() )
        {
            for( int i = 0; i < n - list.size(); i++ )
                list.add( new ArrayList<Double>() );
        }
    }

    private void resizeList(List<Double> list, int n)
    {
        if( n < list.size() )
        {
            for( int i = n; i < list.size(); i++ )
                list.remove( i );
        }
        else if( n > list.size() )
        {
            for( int i = 0; i < n - list.size(); i++ )
                list.add( 0.0 );
        }
    }

    public void sync_to_cycle_model()
    {
        // make sure the inverse map is the right size 
        int n = pCycle_Model.phases.size();
        //        resizeListMap( inverse_index_maps, n );
        resizeListList( transition_rates, n );
        //            inverse_index_maps.resize( n );
        // sync the inverse map to the cell cycle model by 
        // querying the phase_links 
        //            transition_rates.resize( n );
        //        transition_rates.add( new ArrayList<Double>() );
        // also make sure the transition_rates[] are the right size 

        for( int i = 0; i < pCycle_Model.phase_links.size(); i++ )
        {
            //                inverse_index_maps[i].clear(); 
            //            inverse_index_maps.get( i ).clear();
            int size = pCycle_Model.phase_links.get( i ).size();
            for( int j = 0; j < size; j++ )
            {
                //                    inverse_index_maps[i][pCycle_Model.phase_links[i][j].end_phase_index] = j;
                //                    transition_rates[i].resize( pCycle_Model.phase_links[i].size() );
                //j - index of phase_link
                //                int end_index = pCycle_Model.phase_links.get( i ).get( j ).end_phase_index;
                //                inverse_index_maps.get( i ).put( end_index, j );
                List<Double> rates = transition_rates.get( i );
                resizeList( rates, size );
                //                transition_rates.add( i, new ArrayList<Double>( pCycle_Model.phase_links.get( i ).size() ) );
            }
        }
    }

    //        double& Cycle_Data::transition_rate( int start_phase_index , int end_phase_index )
    double transition_rate(int start_phase_index, int end_phase_index)
    {
        //            return transition_rates[ start_phase_index ][ inverse_index_maps[start_phase_index][end_phase_index] ];
        int index = pCycle_Model.inverse_index_maps.get( start_phase_index ).get( end_phase_index );
        return transition_rates.get( start_phase_index ).get( index );

    }

    public void setTransitionRate(int start_phase_index, int end_phase_index, double rate)
    {
        int index = pCycle_Model.inverse_index_maps.get( start_phase_index ).get( end_phase_index );
        transition_rates.get( start_phase_index ).set( index, rate );
    }

    //        double& Cycle_Data::exit_rate(int phase_index )
    double exit_rate(int phase_index)
    {
        //            return transition_rates[phase_index][0];
        return transition_rates.get( phase_index ).get( 0 );
    }
    public void setExitRate(int phase_exit, double rate)
    {
        transition_rates.get( phase_exit ).set( 0, rate );
    }

    public Phase current_phase()
    {
        return pCycle_Model.phases.get( current_phase_index );
    }

    public CycleData clone(CycleModel model)
    {
        CycleData result = new CycleData( model );
        result.time_units = this.time_units;
        result.current_phase_index = this.current_phase_index;
        result.elapsed_time_in_phase = this.elapsed_time_in_phase;
        result.transition_rates = new ArrayList<>();
        for( int i = 0; i < transition_rates.size(); i++ )
        {
            result.transition_rates.add( new ArrayList<Double>( transition_rates.get( i ) ) );
        }
        return result;
    }
}