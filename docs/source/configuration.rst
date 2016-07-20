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
   the mutations are grouped using the column indicated in the `classifier`_.
complement
   same as full but collapsing complementaries
bysample
   the mutations are grouped using the 'SAMPLE' column
file
   The signature is loaded from a file. The optional parameters must be specified.

.. _classifier:

Classifier
   Identifies which column is used as signature identifier. If the column does not exist
   then the classifier value is used as ID.

----

Add inforamtion about the signature and the methods:

From file needs more parameters

Configuration file template
---------------------------
.. literalinclude:: ../../oncodrivefml/oncodrivefml.conf.template
   :language: text
   :linenos:

