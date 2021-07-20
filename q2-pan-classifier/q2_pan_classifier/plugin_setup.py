#   Copyright 2021 Evan Bolyen
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
import importlib

from qiime2.plugin import Str
from qiime2.plugin import Plugin, SemanticType, model

import q2_pan_classifier.actions as actions

from q2_pan_classifier.format_types import MyString, MyStringFormat, MyStringDirFormat

from q2_types.feature_table import FeatureTable, Frequency
# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.
plugin = Plugin("pan_classifier", version="0.0.1.dev",
                website="https://github.com/ebolyen/q2-reveal")

plugin.register_semantic_types(MyString)
plugin.register_semantic_type_to_format(MyString, MyStringDirFormat )
plugin.register_formats(MyStringFormat, MyStringDirFormat)


plugin.methods.register_function(
    function=actions.create_classifier.create_classifier,
    inputs=None,
    outputs=[('cool', MyString)],
    parameters={},
    input_descriptions=None,
    parameter_descriptions=None,
    name='Pan Classifier',
    description="Creates a classifier and classifies them using machine learning"
)

importlib.import_module("q2_pan_classifier.transformers")
