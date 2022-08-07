```css
#toggle {
    position: relative;
    display: inline-block;
    width: 54px;
    height: var(--toggle-height);
    overflow: hidden;
    vertical-align: top;
    overflow: hidden;
    touch-action: none;
}

#toggle input {
    opacity: 0;
    width: 0;
    height: 0;
}

#toggle span {
    position: absolute;
    display: block;
    cursor: pointer;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background-color: #ccc;
    border-radius: var(--toggle-height);
}

#toggle span::before {
    position: absolute;
    content: "";
    height: var(--toggle-slider-radius);
    width: var(--toggle-slider-radius);
    left: var(--toggle-slider-margin);
    bottom: var(--toggle-slider-margin);
    background-color: white;
    border-radius: 50%;
}

#toggle input:checked + span {
    background-color: #48f321;
}

#toggle input:checked + span:before {
    left: 26px;
}
```