package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.SignalBehavior;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class CancerBiorobotsVisualizer extends AgentVisualizer
    {
        @Override
        public Color findBorderColor(Cell cell)
        {
            return Color.black;
        }

        @Override
        public Color findColor(Cell cell)
        {
            Color c = Color.white;
            double damage = SignalBehavior.getSingleSignal( cell, "damage");
            double max_damage = 1.0 * SignalBehavior.getSingleSignal(cell, "custom:damage_rate")
                    / ( 1e-16 + SignalBehavior.getSingleSignal( cell, "custom:repair_rate" ) );

            CellDefinition pCD_cargo = CellDefinition.getCellDefinition( "cargo cell" );
            CellDefinition pCD_cancer = CellDefinition.getCellDefinition( "cancer cell" );
            CellDefinition pCD_worker = CellDefinition.getCellDefinition( "worker cell" );

            //   cargo cell 
            if( cell.type == pCD_cargo.type )
            {
                return Color.BLUE;
                //                output[0] = "blue";
                //                output[1] = "blue";
                //                output[2] = "blue";
                //                output[3] = "none"; // no nuclear outline color 
                //                return output;
            }
                    //  
                    // worker cell 
                    if( cell.type == pCD_worker.type )
                    {
                        return Color.red;
                    }
                    //      output[0] = "red";
                    //      output[1] = "red";
                    //      output[2] = "red"; 
                    //      output[3] = "none"; // no nuclear outline color 
                    //      return output;
                    //  }
                    //  
                    // apoptotic tumor - cyan 
                    if( SignalBehavior.getSingleSignal( cell, "apoptotic" ) > 0.5 ) // Apoptotic - cyan
                    {
                        //                          output[0] = "cyan";
                        //                          output[2] = "darkcyan"; 
                        return Color.cyan;
                    }
                    //  
                    // Necrotic tumor - Brown
                    if( SignalBehavior.getSingleSignal( cell, "necrotic" ) > 0.5 )
                    {
//                        output[0] = "rgb(250,138,38)";
//                        output[2] = "rgb(139,69,19)";
                        return new Color(250,138,18);
                    }
                    //  
                      // live tumor -- shade by level of damage 
                      // if live: color by damage 
                      if( SignalBehavior.getSingleSignal( cell, "dead") < 0.5 )
                      {
                          int damage_int = (int) Math.round( damage * 255.0 / max_damage ); 
                          return new Color( damage_int, 255 - damage_int, damage_int );
                      }
            return c;
        }
    }