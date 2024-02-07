package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class HypothesisRuleset
{

    Map<String, HypothesisRule> rules_map;
    String cell_type;
    CellDefinition pCellDefinition;
    List<HypothesisRule> rules;

    HypothesisRuleset()
    {
        cell_type = "none";
        pCellDefinition = null;
        rules = new ArrayList<>();
        rules_map.clear();
    }

    String display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Behavioral rules for cell type " + cell_type + ":" );
        sb.append( "===================================================" );
        for( HypothesisRule rule : rules )// int i = 0; i < rules.size(); i++ )
        {
            sb.append( rule.reduced_display() );
        }
        sb.append( "\n" );
        return sb.toString();
    }

    String detailed_display()
    {
        StringBuilder sb = new StringBuilder();
        sb.append( "Behavioral rules for cell type " + cell_type + ":" );
        sb.append( "===================================================" );
        for( HypothesisRule rule : rules )// int i = 0; i < rules.size(); i++ )
        {
            sb.append( rule.detailed_display() );
        }
        sb.append( "\n" );
        return sb.toString();
    }


    void sync_to_CellDefinition(CellDefinition pCD)
    {
        pCellDefinition = pCD;
        cell_type = pCD.name;

        for( HypothesisRule rule : rules )// int i = 0; i < rules.size(); i++ )
        {
            rule.sync_to_CellDefinition( pCD );
        }
    }

    HypothesisRule add_behavior( String behavior , double min_behavior, double max_behavior ) throws Exception
    {
        // check: is this a valid signal? (is it in the dictionary?)
        if( SignalBehavior.findBehaviorIndex( behavior ) < 0 )
        {
            throw new Exception( "Warning! Attempted to add behavior " + behavior + " which is not in the dictionary. " + "Either fix your model or add the missing behavior to the simulation."); 

        }

        // first, check. Is there already a ruleset? 
        HypothesisRule rule = rules_map.get( behavior );

            // if not, add it 
//        if( search == rules_map.end() )
//        {
        if (rule == null)
        {
            HypothesisRule hr = new HypothesisRule(); 

            hr.behavior = behavior; 

            hr.sync_to_CellDefinition( pCellDefinition ); 

            hr.min_value = min_behavior; 
            hr.max_value = max_behavior; 

            rules.add( hr );
            //            HypothesisRule pHR = ( rules.back() );
            rules_map.put( behavior, hr );

            return hr;
        }

            // otherwise, edit it 
            //        HypothesisRule pHR = search.second; 

        /*
            // March 28 2023 fix  : let's not overwrite eixsting values
        pHR.min_value = min_behavior; 
        pHR.max_value = max_behavior; 
        */

            return rule;
    }

    HypothesisRule add_behavior(String behavior) throws Exception
    { 
        double min_behavior = 0.1; 
        double max_behavior = 1.0; 
        return add_behavior( behavior, min_behavior, max_behavior );
    }

    void sync_to_CellDefinition(String cell_name)
    {
        sync_to_CellDefinition( CellDefinition.getCellDefinition( cell_name ) );
    }

    HypothesisRule find_behavior(String name)
    {
        return rules_map.get( name );
        //        auto search = rules_map.find( name );
        //        if( search == rules_map.end() )
        //        {
        //            // System.out.println( "Warning! Ruleset does not contain " + name + std::endl; 
        //            // System.out.println( "         Returning NULL." + std::endl; 
        //            return null;
        //        }
        //
        //        return search.second;
    }

    HypothesisRule get(String name)
    {
        return find_behavior( name );
    }

    void apply(Cell pCell) throws Exception
    {
        for( HypothesisRule rule : rules )
        {
            rule.apply( pCell );
        }
    }
}
