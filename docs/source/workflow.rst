Workflow
========

- The elements and the mutations files are loaded (:mod:`oncodrivefml.load`).

    - Load the :ref:`elements <elements dict>`.

    - Generate a tree with the intervals of each element.

    - Load the :ref:`mutations <mutations dict>` and group them by element (using the tree).

- The :ref:`signature <signature dict>` is computed (:mod:`oncodrivefml.signature`).

- Create executors by element (:class:`oncodrivefml.executors.element.ElementExecutor`).

- Launch the executors in parallel using :class:`multiprocessing.pool.Pool`.

- Do a multiple test correction (:mod:`oncodrivefml.mtc`).

- Generate the output (:mod:`oncodrivefml.store`).
