Command Line
============

Below is a list and description of the available command-line tools in ATK. All following commands can be run directly from the command line in any environment in which the package is installed.

   .. _ATKgui:

   **ATKgui**
      opens the ATK gui as if :func:`openGUI() <AstroToolkit.Gui.openGUI>` was called
   
   |

   .. _ATKquery:

   **ATKquery**
      runs an ATK query as if :func:`query() <AstroToolkit.Tools.query>` was called. The standard :ref:`data structure <Data Structures>` methods are available as "jobs" for any returned data

      **Usage:**
         .. code-block:: python
            
            ATKquery kind survey source/pos radius/size

      Providing a **radius**/**size** is optional, with the default given by the relevant :ref:`config key <Config Keys>`. See :func:`query() <AstroToolkit.Tools.query>` for a description of these parameters in each kind of query.
      
      **Note:** 
         When using the savedata, showplot and saveplot jobs on returned data, an additional argument maybe be entered which will provide an override to the file name. See the savedata method in :ref:`Data Structures` for more information.

         This tool is provided for quick command-line queries, and is therefore more limited in functionality compared to queries performed through either the :ref:`ATK GUI <Gui>` or scripts.

         To exit the ATK command-line interface, enter 'exit' when choosing a job.

   |
    
   .. _ATKsearch:

   **ATKsearch**
      searches for a target in Vizier/SIMBAD as if :func:`search <AstroToolkit.Tools.search>` was called.

      **Usage:**
         .. code-block:: python

            ATKsearch kind source/pos radius

      Providing a radius is optional, with the default given by the :ref:`search_radius <cfg_search_radius>` config key. See :func:`search <AstroToolkit.Tools.search>` for a description of the available parameters.

   |

   .. _ATKread:

   **ATKread**
      reads a local ATK data file as if :func:`readdata <AstroToolkit.Tools.readdata>` was called.

      **Usage:**
         .. code-block:: python

            ATkread fname

      where **fname** is the name of the ATK file to be read.

      **Note:** 
         When using the savedata, showplot and saveplot jobs on returned data, an additional argument maybe be entered which will provide an override to the file name. See the savedata method in :ref:`Data Structures` for more information.

         To exit the ATK command-line interface, enter 'exit' when asked which job to perform.

   |

   .. _ATKshowconfig:

   **ATKshowconfig**
      prints the ATK config to stdout
    
   |

   .. _ATKeditconfig:

   **ATKeditconfig**
      sets the value of a given key in the ATK config

      **Usage:**
         .. code-block:: console
            
            ATKeditconfig <key> <value>

      where **key** is the name of the config key and **value** is its desired value

   |
    
   .. _ATKopenconfig:

   **ATKopenconfig**
      opens the ATK config in the default text editor

   |

   .. _ATKresetconfig:

   **ATKresetconfig**
      resets the ATK config to default values

   |

   .. _ATKshowaliases:
    
   **ATKshowaliases**
      prints all catalogue aliases to stdout

   |

   .. _ATKaddalias:

   **ATKaddalias**
      adds a catalogue alias
      
      **Usage:**
         .. code-block:: console
            
            ATKaddalias <name> <id>

      where **name** is the name of the alias (e.g. allwise) and **id** is its Vizier catalogue ID (e.g. II/328/allwise) 

   |

   .. _ATKdelalias:

   **ATKdelalias**
      removes an existing catalogue alias
      
      **Usage:**
         .. code-block:: console

            ATKdelalias <name>

      where **name** is the name of the alias

   |

   .. _ATKopenaliases:

   **ATKopenaliases**
      opens the ATK catalogue aliases file in the default text editor

   |

   .. _ATKresetaliases:

   **ATKresetaliases**
      resets the alias file to its default state (i.e. deletes all added aliases)

   |

   .. _ATKconv:

   **ATKconv** 
      converts HMS±DMS coordinates to degrees and vice versa

      **Usage:**
         .. code-block:: python

            ATKconv position

      where **position** is either two arguments (ra and dec in degrees), or a single position in HMS±DMS format. See :func:`hms2deg <AstroToolkit.Tools.hms2deg>` and :func:`deg2hms <AstroToolkit.Tools.deg2hms>` for reference.
