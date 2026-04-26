"""
##########
The Config
##########
ATK uses a config file to manage many default behaviours. **For a list and description of all available config keys, see the** :doc:`next tutorial <config_keys>`.

Opening the Config
==================
The config can be opened in the default text editor, either from the commandline:

.. code-block:: console

    $ ATKconfig open

or from inside a script:

.. code-block:: python

    from ATK.Config import CONFIG

    CONFIG.open()

|

Viewing the Config
==================
To print the config to the terminal, use:

.. code-block:: console

    $ ATKconfig show

or:
"""

from ATK.Config import CONFIG

CONFIG.show()

# sphinx_gallery_start_ignore
from ATK.configuration.defaults.base_defaults import CONFIG_KEY_DEFS, DEFAULTS

with open("../../auto_tutorials/configuration/keys.rst", "w") as f:
    for section in DEFAULTS:
        f.write(f".. rubric:: {section}\n\n")
        for key in DEFAULTS[section]:
            desc = CONFIG_KEY_DEFS[section].get(key)
            f.write(f"- **{key}**\n")
            f.write(f"      {desc}")
            f.write(f" Default is {DEFAULTS[section][key]}.\n")
        f.write("\n|\n\n")
# sphinx_gallery_end_ignore

# %%
#
# |
#
# Editing the config
# ==================
# The value of a key in a given section of the config can be set as follows:
#
# .. code-block:: console
#
#    $ ATKconfig set <section> <key> <value>
#
# .. code-block:: python
#
#    from ATK.Config import CONFIG
#
#    CONFIG[section][key] = [value]
#
# |
#
# Resetting the Config
# ====================
# The config can also be reset to its default state:
#
# .. code-block:: console
#
#    $ ATKconfig reset
#
# .. code-block:: python
#
#    from ATK.Config import CONFIG
#
#    CONFIG.reset()

# %%
#
# |
# |
# |
#
# .. rubric:: Download this Tutorial
