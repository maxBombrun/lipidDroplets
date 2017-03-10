# lipidDroplets

Detect individual lipid droplets, usually clustered tightly.


Works with screening from BODIPY and HOESCHT channel.

Developped in Python.

## What it does

The program takes a screening in BODIPY and HOESCHT channels as input and provide images of the segmentation of individual lipid droplets and associated nuclei.
It also extracts several features such as the size distribution of the individual droplets per cell.

The aim is to improve models for High-Content/High-Throughput Microscopy Analysis of Subject-Specific Adipogenesis.



## Usage

Modify the settings.py file to define the correct path, mostly input data.
In the main.py, you can chose the modules you want to apply.

On the first run, apply them sequentially.

This algorithm needs CellProfiler to define the relationship between nuclei and associated lipid droplets, i.e., approximated the cells.

http://cellprofiler.org/releases/



## Credits

This project was funded by Science for Life Laboratory, the Swedish research council (grant 2012-4968 to Carolina Wählby and K2012-55X-010334-46-5-to Peter Arner), the Swedish strategic research program eSSENCE (to Carolina Wählby) and CIMED (to Peter Arner).

## License 

    Copyright (c) 2016-2017, Maxime Bombrun
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification,
    are permitted provided that the following conditions are met:

      Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

      Redistributions in binary form must reproduce the above copyright notice, this
      list of conditions and the following disclaimer in the documentation and/or
      other materials provided with the distribution.

      Neither the name of the {organization} nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
