"""
The Config
==========

Opening the Config
------------------
.. code-block:: console

    $ ATKconfig open

or:

.. code-block:: python

    from ATK.Config import CONFIG

    Config.open()

|

Viewing the Config
-------------------

.. code-block:: console

    $ ATKconfig open

or:
"""

from ATK.Config import CONFIG

CONFIG.show()

# %%
#
# |
#
# Editing the config
# ------------------
#
# .. code-block:: console
#
#    $ ATKconfig set <section> <key> <value>
#
# or:
#
# .. code-block:: python
#
#    from ATK.Config import CONFIG
#    CONFIG.<section>.<key> = <value>

# %%
# |
# |
# |
#
# Download this Tutorial
# ----------------------
