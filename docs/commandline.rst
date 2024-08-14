Command-Line
============

Here, command-line integration will be outlined. All following commands can be run directly from the command-line.

   **ATKgui**
      opens the ATK gui as if Gui.openGUI was called
   
   |

   **ATKtsguide**
      outputs the PyAOV installation instructions for your platform as if Setup.tsguide was called

   |

   **ATKtsbuild**
      builds PyAOV tools for your platform as if Setup.tsbuild was called

   |

   **ATKquery**
      runs an ATK query as if Tools.query was called, and provides the showdata, savedata, showplot and saveplot methods as "jobs"

      **Arguments:**

         ATKquery kind survey source/pos radius/size

      where providing a radius/size is optional. See 'query' in :ref:`Tools` for a description of these parameters in each kind of query.
      
      **Note:** 
         When using the savedata, showplot and saveplot jobs on returned data, an additional argument maybe be entered which will provide an override to the file name. See the savedata method in :ref:`Data Structures` for more information.

         This tool is provided for quick command-line queries, and is therefore more limited in functionality compared to queries performed through either the ATK GUI or scripts. Some additional configuration may be possible through the :ref:`Config`.

         To exit the ATK command-line interface, enter 'exit' when asked which job to perform.

   |

   **ATKread**
      reads a local ATK file as if Tools.readdata was called, and provides the showdata, savedata, showplot and saveplot methods as "jobs"

      **Arguments:**

         ATkread fname

      where fname is the name of the ATK file to be read.

      **Note:** 
         When using the savedata, showplot and saveplot jobs on returned data, an additional argument maybe be entered which will provide an override to the file name. See the savedata method in :ref:`Data Structures` for more information.

         This tool is provided for quick command-line queries, and is therefore more limited in functionality compared to queries performed through either the ATK GUI or scripts. Some additional configuration may be possible through the :ref:`Config`.

         To exit the ATK command-line interface, enter 'exit' when asked which job to perform.

   |
