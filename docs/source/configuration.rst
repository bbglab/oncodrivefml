.. _project configuration:

Configuration
=============

The method behaviour can be modified through a configuration file.

Specifications
--------------

That configuration file has to fulfill these specifications:

.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template.spec
   :language: text
   :linenos:

Any field without a default value must appear.


Signature
^^^^^^^^^

none
   :ref:`signature dict <signature dict>` is None and all the probabilities for all mutations are the same.
full

complement
   same as full but collapsing complementaries
bysample

file
   The signature is loaded from a file. The optional parameters must be specified.



----

Add inforamtion about the signature and the methods:

From file needs more parameters

Configuration file template
---------------------------
.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :linenos:

